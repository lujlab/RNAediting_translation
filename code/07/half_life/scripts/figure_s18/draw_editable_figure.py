#!/usr/bin/env python3
from pathlib import Path
import json,math,re
import fitz

import argparse
parser=argparse.ArgumentParser(description='Render editable Figure S18 from recomputed statistics.')
parser.add_argument('--data-dir',type=Path,required=True)
parser.add_argument('--output-pdf',type=Path,required=True)
args=parser.parse_args()
OUT=args.data_dir.resolve()
ASSETS=Path(__file__).resolve().parent/'assets'
layout=json.loads((ASSETS/'layout.json').read_text(encoding='utf-8'))
data=json.loads((OUT/'plot_data.json').read_text(encoding='utf-8'))
doc=fitz.open();page=doc.new_page(width=layout['width'],height=layout['height'])
BLACK=(0,0,0);RED=(.894,.341,.337);BLUE=(.298,.471,.659);GREEN=(.329,.635,.294)
fonts={};fontobjects={}
for name,filename in [('regular','Arial_regular.ttf'),('italic','Arial_italic.ttf'),('bold','Arial_bold.ttf'),('math','Arial_math.ttf')]:
    buffer=(ASSETS/filename).read_bytes()
    fonts[name]=buffer;fontobjects[name]=fitz.Font(fontbuffer=buffer)
    page.insert_font(fontname=name,fontbuffer=buffer)


def text(x,y,s,size=7,font='regular',anchor='left',rotate=0):
    s=s.replace('\u00ad','-')
    assert all(fontobjects[font].has_glyph(ord(c)) for c in s),repr(s)
    width=fontobjects[font].text_length(s,fontsize=size)
    if not rotate:
        if anchor=='center':x-=width/2
        if anchor=='right':x-=width
    page.insert_text((x,y),s,fontname=font,fontsize=size,rotate=rotate,color=BLACK)
    return width


def rich(x,y,runs,size=7,anchor='left'):
    width=sum(fontobjects[f].text_length(s,fontsize=size) for s,f in runs)
    if anchor=='center':x-=width/2
    elif anchor=='right':x-=width
    for s,f in runs:x+=text(x,y,s,size,f)


def line(x1,y1,x2,y2,color=BLACK,width=.5,dashes=None):
    page.draw_line((x1,y1),(x2,y2),color=color,width=width,dashes=dashes)


def polygon(points,fill,stroke=None,width=.5):
    shape=page.new_shape();shape.draw_polyline(points)
    shape.finish(color=stroke,fill=fill,width=width,closePath=True)
    shape.commit()


def bounded_polygon(points,top,bottom):
    for bound,above in [(top,True),(bottom,False)]:
        output=[]
        for a,b in zip(points[-1:]+points[:-1],points):
            ia=(a[1]>=bound) if above else (a[1]<=bound)
            ib=(b[1]>=bound) if above else (b[1]<=bound)
            if ia!=ib:
                fraction=(bound-a[1])/(b[1]-a[1])
                output.append((a[0]+fraction*(b[0]-a[0]),bound))
            if ib:output.append(b)
        points=output
    return points


def rect(x1,y1,x2,y2,fill,stroke=None,width=.5):
    page.draw_rect(fitz.Rect(x1,y1,x2,y2),fill=fill,color=stroke,width=width)


def circle(x,y,r,color):
    page.draw_circle((x,y),r,color=None,fill=color)


def diamond(x,y):
    polygon([(x,y-4),(x+4,y),(x,y+4),(x-4,y)],(1,1,1),BLACK,.5)


def number_n(x,y,n):
    rich(x,y,[('n','italic'),(' = {:,}'.format(n),'regular')],7,anchor='center')


def violin_panel(groups,centers,ymap,ticks,axes,label_positions,half=False):
    left,right,top,bottom=axes
    colorlist=[RED,BLUE];spacing=centers[1]-centers[0]
    maximum=max(max(g['density']) for g in groups)
    if not half:line(left,ymap(math.log10(3000)),right,ymap(math.log10(3000)),(.87,.87,.87),.5)
    for g,xc,color in zip(groups,centers,colorlist):
        ys=g['density_y'];density=g['density']
        pts=[(xc+.45*spacing*v/maximum,ymap(y)) for y,v in zip(ys,density)]
        pts += [(xc-.45*spacing*v/maximum,ymap(y)) for y,v in reversed(list(zip(ys,density)))]
        polygon(bounded_polygon(pts,top,bottom),color,BLACK,.5)
        bw=.18*spacing
        line(xc,ymap(g['whisk_lo']),xc,ymap(g['whisk_hi']),BLACK,.5)
        rect(xc-bw/2,ymap(g['q3']),xc+bw/2,ymap(g['q1']),color,BLACK,.5)
        line(xc-bw/2,ymap(g['median']),xc+bw/2,ymap(g['median']),BLACK,.5)
        lower=max(5 if half else 200,g['mean']-g['sd'])
        upper=g['mean']+g['sd']
        y1,y2=ymap(math.log10(lower)),ymap(math.log10(upper))
        line(xc,y1,xc,y2,BLACK,.5)
        for yy in [y1,y2]:line(xc-.045*spacing,yy,xc+.045*spacing,yy,BLACK,.5)
        diamond(xc,ymap(math.log10(g['mean'])))
        value=('{:.1f} ± {:.1f}'.format(g['mean'],g['sd']) if half else '{:,.0f} ± {:,.0f}'.format(g['mean'],g['sd']))
        text(xc,33.95,value,6,anchor='center')
        number_n(xc,138.4,g['n'])
    line(left,top,left,bottom);line(left,bottom,right,bottom)
    for val,yy in ticks:line(left-1.5,yy,left,yy)
    for xc in centers:line(xc,bottom,xc,bottom+2.3)


# A. Transcript lengths
Ay=lambda z:110.71-33.09*(z-math.log10(300))
Aticks=[(x,Ay(math.log10(x))) for x in [300,1000,3000,10000,30000]]
violin_panel(data['A'],[70.6213,124.8758],Ay,Aticks,(38.08,157.417,22.61,120.85),None)

# B. Decay classes
b=data['B'];base,bar_top=121.25,41.56
for xc,width,dec,stable in [(219.77,33.66,b['edited_decay'],b['edited_stable']),
                           (266.51,33.66,b['control_decay'],b['control_stable'])]:
    n=dec+stable;frac=stable/n;boundary=base-(base-bar_top)*frac
    rect(xc-width/2,bar_top,xc+width/2,boundary,RED)
    rect(xc-width/2,boundary,xc+width/2,base,BLUE)
    for count,ycenter in [(dec,(bar_top+boundary)/2),(stable,(boundary+base)/2)]:
        rich(xc,ycenter-1,[('n','italic'),(' = {:,}'.format(count),'regular')],6.5,anchor='center')
        text(xc,ycenter+7,'{:.1f}%'.format(100*count/n),6.5,anchor='center')
line(191.72,28.01,191.72,121.25);line(191.72,121.25,294.56,121.25)
for val in [0,25,50,75,100]:
    yy=base-(base-bar_top)*val/100
    line(190.22,yy,191.72,yy)
for xc in [219.77,266.51]:line(xc,121.25,xc,123.55)
line(214.16,36.38,272.12,36.38)
line(214.16,36.38,214.16,38.55);line(272.12,36.38,272.12,38.55)
rich(243.14,33.95,[('P','italic'),(' = {:.3g}'.format(b['P']),'regular')],6.5,anchor='center')
rect(210.88,14.61,217.79,21.52,RED);rect(245.61,14.61,252.52,21.52,BLUE)

# C. mRNA half-lives
Cy=lambda z:109.82-22.332*(z-1)
Cticks=[(x,Cy(math.log10(x))) for x in [10,100,1000,10000]]
violin_panel(data['C'],[372.0881,427.2286],Cy,Cticks,(339,457.743,22.61,120.85),None,True)

# D. Regression estimates
order=['decay_classIV','edited','z_A_content','z_GC_content','z_log2_UTR5_length',
       'z_log2_UTR3_length','z_log2_CDS_length','z_log2_transcript_length',
       'z_log2_published_expression','decay_classV','decay_classIII']
index={v:i for i,v in enumerate(order)}
dx=lambda z:197.68+65.877*z
line(117.31,171.68,117.31,277.404);line(117.31,277.234,283.31,277.234)
line(dx(0),173.782,dx(0),277.234,BLACK,.35,'[4 4] 0')
for val in [-1,-.5,0,.5,1]:line(dx(val),277.234,dx(val),279.55)
colors={'H3':RED,'H4':GREEN,'H5':BLUE};offset={'H3':1.316,'H4':0,'H5':-1.316}
for row in data['D']:
    if row['variable'] not in index:continue
    assert -1.22<=row['lower']<=row['upper']<=1.30, row
    yy=177.2899+9.401*index[row['variable']]+offset[row['model']]
    color=colors[row['model']]
    line(dx(row['lower']),yy,dx(row['upper']),yy,color,.4)
    circle(dx(row['estimate']),yy,1.05,color)
for nm,y in [('H3',250.97),('H4',259.36),('H5',267.75)]:
    line(214.41,y,219.67,y,colors[nm],.4)

# E. CDS editing level and mRNA half-life
e=data['E'];ex=lambda x:334.2193+(x+18)*13.28195;ey=lambda y:269.0408-(y-4)*13.23365
band=[(ex(x),ey(y)) for x,y in zip(e['fit_x'],e['upper'])]
band += [(ex(x),ey(y)) for x,y in reversed(list(zip(e['fit_x'],e['lower'])))]
assert all(330.8<=x<=457.9 and 171.1<=y<=277 for x,y in band)
polygon(band,(.8835,.8835,.8835))
for x1,y1,x2,y2 in zip(e['fit_x'][:-1],e['fit'][:-1],e['fit_x'][1:],e['fit'][1:]):
    line(ex(x1),ey(y1),ex(x2),ey(y2),(.651,.239,.251),.4)
for x,y in zip(e['x'],e['y']):circle(ex(x),ey(y),1.4,RED)
line(330.9463,171.2263,330.9463,276.7843);line(330.9463,276.7843,457.7426,276.7843)
for y in range(4,12):line(329.1183,ey(y),330.9463,ey(y))
for x in range(-18,-8):line(ex(x),276.7843,ex(x),278.4793)
rich(333.7002,160.1602,[("Spearman's ",'regular'),('ρ','italic'),
    (' = {:.3f},  '.format(e['rho']).replace('-','−'),'regular'),('P','italic'),
    (' = {:.3f}'.format(e['P']),'regular')],7)
rich(333.7002,168.55,[('n','italic'),(' = {}'.format(e['n']),'regular')],7)

static_count=0
for label in layout['labels']:
    tx,ty=label['x'],label['y']
    for part,face in label['runs']:
        advance=text(tx,ty,part,label['size'],face,rotate=90 if label['vertical'] else 0)
        if label['vertical']:ty-=advance
        else:tx+=advance
    static_count+=1

nrow=next(r for r in data['D'] if r['model']=='H3')
rich(119.10,160.16,[('N','italic'),(' = {:,} (Edited = {:,}; control = {:,})'.format(
     nrow['N'],nrow['N_edited'],nrow['N_unedited']),'regular')],7)

doc.set_metadata(dict(title='Supplementary Figure 18 — r6.12 GTF CDS transcripts',
                     subject='Native editable vector figure; FlyBase r6.12 GTF / dm6; longest CDS among site-containing transcripts',
                     author='RNA editing analysis',creator='Native vector paths and editable embedded Arial text'))
combined=b'\n'.join(doc.xref_stream(x) for x in page.get_contents())
stream=page.get_contents()[0];doc.update_stream(stream,combined)
page.set_contents(stream)
path=args.output_pdf.resolve()
path.parent.mkdir(parents=True,exist_ok=True)
doc.save(str(path),garbage=4,deflate=True)
doc.close()

check=fitz.open(str(path));p=check[0]
streams=b'\n'.join(check.xref_stream(x) or b'' for x in p.get_contents())
objects='\n'.join(check.xref_object(i) for i in range(1,check.xref_length()))
report=dict(page_width=p.rect.width,page_height=p.rect.height,pages=len(check),
    page_content_streams=len(p.get_contents()),
    images=len(p.get_images(full=True)),form_xobjects=len(p.get_xobjects()),
    clipping_operators=len(re.findall(rb'(?<!\S)W\*?(?!\S)',streams)),
    masks=bool(re.search(r'/(?:SMask|Mask)\b',objects)),transparency_groups='/Group' in objects,
    optional_layers='/OCProperties' in objects,editable_text_characters=len(p.get_text()),
    static_label_spans=static_count,point_count=e['n'])
assert report['images']==report['form_xobjects']==report['clipping_operators']==0
assert not report['masks'] and not report['transparency_groups'] and not report['optional_layers']
assert report['editable_text_characters']>500
assert not re.search(r'n\s*=\s*147',p.get_text()),p.get_text()
(OUT/'PDF_editability_audit.json').write_text(json.dumps(report,indent=2))
p.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False).save(str(OUT/'Supplementary_Figure_18_r6_12_transcripts_preview.png'))
print(path)
print(json.dumps(report,indent=2))
