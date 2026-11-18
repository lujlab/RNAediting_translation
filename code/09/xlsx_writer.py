import zipfile, xml.etree.ElementTree as ET
from xml.sax.saxutils import escape


def write_workbook(path, sheets):
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"

    def col(n):
        s = ""
        while n:
            n, r = divmod(n - 1, 26)
            s = chr(65 + r) + s
        return s

    def sheetxml(head, data):
        allrows = [head] + data
        x = [
            f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><cols>'
        ]
        for i in range(len(head)):
            x.append(
                f'<col min="{i+1}" max="{i+1}" width="{85 if len(head)==2 and i==1 else 30}" customWidth="1"/>'
            )
        x.append("</cols><sheetData>")
        for ri, row in enumerate(allrows, 1):
            x.append(f'<row r="{ri}" ht="{60 if ri==1 else 44}" customHeight="1">')
            for ci, value in enumerate(row, 1):
                ref = f"{col(ci)}{ri}"
                style = 1 if ri == 1 else 0
                if isinstance(value, (int, float)):
                    x.append(f'<c r="{ref}" s="{style}"><v>{value}</v></c>')
                else:
                    x.append(
                        f'<c r="{ref}" s="{style}" t="inlineStr"><is><t xml:space="preserve">{escape(str(value))}</t></is></c>'
                    )
            x.append("</row>")
        x.append(
            f'</sheetData><autoFilter ref="A1:{col(len(head))}{len(allrows)}"/></worksheet>'
        )
        return "".join(x)

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        types = [
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        ]
        for i in range(1, len(sheets) + 1):
            types.append(
                f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            )
        z.writestr("[Content_Types].xml", "".join(types) + "</Types>")
        z.writestr(
            "_rels/.rels",
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        )
        z.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
            + "".join(
                f'<sheet name="{name}" sheetId="{i}" r:id="rId{i}"/>'
                for i, (name, _, _) in enumerate(sheets, 1)
            )
            + "</sheets></workbook>",
        )
        rels = [
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        ]
        for i, (name, head, data) in enumerate(sheets, 1):
            rels.append(
                f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
            )
            z.writestr(f"xl/worksheets/sheet{i}.xml", sheetxml(head, data))
        rels.append(
            '<Relationship Id="styles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
        )
        z.writestr("xl/_rels/workbook.xml.rels", "".join(rels))
        z.writestr(
            "xl/styles.xml",
            f'<styleSheet xmlns="{ns}"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>',
        )
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        for name in z.namelist():
            ET.fromstring(z.read(name))
        actual = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        assert len(actual.findall(f".//{{{ns}}}row")) == len(sheets[0][2]) + 1
