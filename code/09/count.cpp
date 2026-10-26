#include <bits/stdc++.h>
using namespace std;
struct Block {
    int lo, hi, t0;
};
struct Candidate {
    int mask;
    string ref;
};
struct Tx {
    string name, chr;
    bool minus;
    int cs, ce, chrid;
    vector<Block> b;
    unordered_map<int, Candidate> cand;
};
struct Count {
    double dp = 0, alt = 0, allA = 0;
    double A[3] = {0, 0, 0}, G[3] = {0, 0, 0};
};
struct Key {
    uint32_t chrom, p0, p1, p2;
    bool operator==(const Key &o) const {
        return chrom == o.chrom && p0 == o.p0 && p1 == o.p1 && p2 == o.p2;
    }
};
struct KH {
    size_t operator()(const Key &k) const {
        size_t h = k.chrom;
        for (auto x : {k.p0, k.p1, k.p2})
            h ^= size_t(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
struct CInfo {
    double dp = 0, alt = 0;
    bool pass = false;
    string besttx;
};
vector<string> split(const string &s, char c = '\t') {
    vector<string> a;
    size_t p = 0, q;
    while ((q = s.find(c, p)) != string::npos) {
        a.push_back(s.substr(p, q - p));
        p = q + 1;
    }
    a.push_back(s.substr(p));
    return a;
}
int gpos(const Tx &t, int p) {
    for (const auto &b : t.b)
        if (p >= b.t0 && p < b.t0 + b.hi - b.lo + 1)
            return t.minus ? b.hi - (p - b.t0) : b.lo + (p - b.t0);
    return -1;
}
Key key(const Tx &t, int p) {
    int a = gpos(t, p), b = gpos(t, p + 1), c = gpos(t, p + 2);
    if (a < 0 || b < 0 || c < 0)
        throw runtime_error("unmapped codon " + t.name + ":" + to_string(p));
    return Key{uint32_t(t.chrid * 2 + t.minus), uint32_t(a), uint32_t(b), uint32_t(c)};
}
int main(int argc, char **argv) {
    try {
        if (argc != 5)
            throw runtime_error("usage: count OUT SPECIES DATASET LIBRARY");
        string out = argv[1], sp = argv[2], ds = argv[3], lib = argv[4], s;
        vector<Tx> tx;
        unordered_map<string, int> ids, chromids;
        vector<string> chromnames;
        ifstream meta(out + "/" + sp + ".transcripts.tsv");
        while (getline(meta, s)) {
            auto f = split(s);
            Tx t;
            t.name = f[0];
            t.chr = f[1];
            t.minus = f[2] == "-";
            t.cs = stoi(f[3]);
            t.ce = stoi(f[4]);
            if (!chromids.count(t.chr)) {
                chromids[t.chr] = chromnames.size();
                chromnames.push_back(t.chr);
            }
            t.chrid = chromids[t.chr];
            int off = 1;
            for (auto x : split(f[5], ',')) {
                auto z = split(x, '-');
                int a = stoi(z[0]), b = stoi(z[1]);
                t.b.push_back({a, b, off});
                off += b - a + 1;
            }
            ids[t.name] = tx.size();
            tx.push_back(move(t));
        }
        if (tx.empty())
            throw runtime_error("empty transcript metadata");
        ifstream cand(out + "/" + sp + ".candidates.tsv");
        getline(cand, s);
        while (getline(cand, s)) {
            auto f = split(s);
            tx.at(ids.at(f[0])).cand[stoi(f[1])] = {stoi(f[2]), f[3]};
        }
        vector<string> paths;
        ifstream manifest(out + "/inputs.tsv");
        getline(manifest, s);
        while (getline(manifest, s)) {
            if (!s.empty() && s.back() == '\r')
                s.pop_back();
            auto f = split(s);
            if (f[0] == sp && f[1] == ds && f[2] == lib)
                paths.push_back(f[4]);
        }
        if (paths.empty())
            throw runtime_error("no inputs");
        ofstream mismatches(out + "/" + sp + "." + ds + "." + lib +
                            ".coordinate_checks.tsv");
        mismatches << "input_file\ttranscript\ttranscript_codon_start\trecorded_"
                      "chromosome\trecorded_strand\trecorded_base1\treference_"
                      "chromosome\treference_strand\treference_base1\n";
        unordered_map<uint64_t, Count> counts;
        counts.reserve(100000);
        map<string, uint64_t> qc;
        for (auto path : paths) {
            ifstream in(path);
            if (!in)
                throw runtime_error("cannot open " + path);
            string prev;
            int tid = -1;
            uint64_t lines = 0;
            while (getline(in, s)) {
                ++lines;
                array<string_view, 7> f;
                size_t p = 0;
                bool bad = false;
                for (int j = 0; j < 7; j++) {
                    auto q = s.find('\t', p);
                    if (j < 6 && q == string::npos) {
                        bad = true;
                        break;
                    }
                    f[j] = string_view(s.data() + p,
                                       (q == string::npos ? s.size() : q) - p);
                    p = q + 1;
                }
                if (bad) {
                    ++qc["bad_columns"];
                    continue;
                }
                if (f[0] != prev) {
                    prev = string(f[0]);
                    auto it = ids.find(prev);
                    tid = it == ids.end() ? -1 : it->second;
                }
                if (tid < 0) {
                    ++qc["unknown_transcript"];
                    continue;
                }
                int pos = atoi(f[1].data());
                double n = strtod(f[3].data(), nullptr);
                if (!(n > 0) || !isfinite(n)) {
                    ++qc["nonpositive_or_invalid_count"];
                    continue;
                }
                const Tx &t = tx[tid];
                if (!t.cand.count(pos))
                    continue;
                if (pos < t.cs || pos + 2 > t.ce) {
                    ++qc["outside_complete_CDS"];
                    continue;
                }
                if ((pos - t.cs) % 3) {
                    ++qc["out_of_frame"];
                    continue;
                }
                if (f[2].size() != 3) {
                    ++qc["bad_codon_length"];
                    continue;
                }
                // Verify recorded genomic first base against the exon-aware mapping.
                if (gpos(t, pos) != atoi(f[6].data()) || f[4] != t.chr ||
                    f[5] != (t.minus ? "-" : "+")) {
                    ++qc["coordinate_mismatch"];
                    mismatches << path << '\t' << t.name << '\t' << pos << '\t' << f[4]
                               << '\t' << f[5] << '\t' << f[6] << '\t' << t.chr << '\t'
                               << (t.minus ? "-" : "+") << '\t' << gpos(t, pos) << '\n';
                    continue;
                }
                auto &c = counts[(uint64_t(tid) << 32) | uint32_t(pos)];
                c.dp += n;
                auto ci = t.cand.find(pos);
                if (ci != t.cand.end()) {
                    bool alt = false;
                    for (int j = 0; j < 3; j++)
                        if ((ci->second.mask & (1 << j)) && f[2][j] == 'G')
                            alt = true;
                    if (alt)
                        c.alt += n;
                    bool allA = true;
                    for (int j = 0; j < 3; j++) {
                        if (f[2][j] == 'A')
                            c.A[j] += n;
                        if (f[2][j] == 'G')
                            c.G[j] += n;
                        if ((ci->second.mask & (1 << j)) && f[2][j] != 'A')
                            allA = false;
                    }
                    if (allA)
                        c.allA += n;
                }
                ++qc["included_codon_state_rows"];
            }
            qc["input_rows"] += lines;
            cerr << sp << " " << ds << " " << lib << " read " << path
                 << " rows=" << lines << " pooled_tx_codons=" << counts.size() << endl;
        }
        string prefix = out + "/" + sp + "." + ds + "." + lib;
        ofstream detail(prefix + ".transcript_candidate_counts.tsv");
        detail << "transcript\tcodon_start_"
                  "tx\tchromosome\tstrand\tbase1\tbase2\tbase3\tref_codon\tedit_"
                  "mask\ttotal\talt\tallA\tA1\tG1\tA2\tG2\tA3\tG3\n"
               << setprecision(17);
        for (int tid = 0; tid < (int)tx.size(); tid++) {
            auto &t = tx[tid];
            for (auto &ce : t.cand) {
                int p = ce.first;
                auto k = key(t, p);
                auto v = counts[(uint64_t(tid) << 32) | uint32_t(p)];
                detail << t.name << '\t' << p << '\t' << t.chr << '\t'
                       << (t.minus ? "-" : "+") << '\t' << k.p0 << '\t' << k.p1 << '\t'
                       << k.p2 << '\t' << ce.second.ref << '\t' << ce.second.mask
                       << '\t' << v.dp << '\t' << v.alt << '\t' << v.allA;
                for (int j = 0; j < 3; j++)
                    detail << '\t' << v.A[j] << '\t' << v.G[j];
                detail << '\n';
            }
        }
        cerr << "DONE " << prefix << " candidate transcript codons=" << counts.size()
             << endl;
    } catch (exception &e) {
        cerr << "ERROR " << e.what() << endl;
        return 1;
    }
    return 0;
}
