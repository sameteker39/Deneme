const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle, PageBreak
} = require('docx');

const FONT = 'Cambria';
const SIZE = 24;        // 12pt (half-points)
const TOTAL_W = 9360;   // usable width in DXA for A4 with 1" margins

// ---- inline markdown -> TextRun[] -------------------------------------------
function inline(text, base = {}) {
  const runs = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)/g;
  let last = 0, m;
  const push = (t, opts) => { if (t) runs.push(new TextRun({ text: t, font: FONT, size: SIZE, ...base, ...opts })); };
  while ((m = re.exec(text)) !== null) {
    push(text.slice(last, m.index));
    const tok = m[0];
    if (tok.startsWith('**')) push(tok.slice(2, -2), { bold: true });
    else if (tok.startsWith('`')) push(tok.slice(1, -1), { font: 'Consolas', size: SIZE - 2 });
    else push(tok.slice(1, -1), { italics: true });
    last = m.index + tok.length;
  }
  push(text.slice(last));
  return runs.length ? runs : [new TextRun({ text: '', font: FONT, size: SIZE, ...base })];
}

const para = (text, o = {}) => new Paragraph({
  children: inline(text, o.run || {}),
  spacing: { line: o.line ?? 360, after: o.after ?? 180 },
  alignment: o.align,
  indent: o.indent,
  heading: o.heading,
  numbering: o.numbering,
  bullet: o.bullet,
  border: o.border,
});

function splitRow(line) {
  return line.trim().replace(/^\|/, '').replace(/\|$/, '').split('|').map(s => s.trim());
}

function makeTable(rows) {
  const header = splitRow(rows[0]);
  const body = rows.slice(2).map(splitRow);
  const n = header.length;
  const colW = Math.floor(TOTAL_W / n);
  const widths = Array(n).fill(colW);
  widths[n - 1] = TOTAL_W - colW * (n - 1);

  const cell = (txt, isHead) => new TableCell({
    width: { size: widths[0], type: WidthType.DXA },
    shading: isHead ? { type: ShadingType.CLEAR, fill: 'EDEDED' } : undefined,
    margins: { top: 60, bottom: 60, left: 90, right: 90 },
    children: [new Paragraph({
      children: inline(txt, isHead ? { bold: true } : {}),
      spacing: { line: 260, after: 0 },
    })],
  });

  const mkRow = (cells, isHead) => new TableRow({
    tableHeader: isHead,
    children: cells.map((c, i) => {
      const tc = cell(c, isHead);
      tc.root[1].root.width = { size: widths[i], type: WidthType.DXA };
      return tc;
    }),
  });

  return new Table({
    columnWidths: widths,
    width: { size: TOTAL_W, type: WidthType.DXA },
    rows: [mkRow(header, true), ...body.map(r => mkRow(r, false))],
  });
}

// ---- block parser -----------------------------------------------------------
function parse(md) {
  const out = [];
  const lines = md.split('\n');
  let i = 0;
  while (i < lines.length) {
    const line = lines[i];

    if (!line.trim()) { i++; continue; }

    // table
    if (/^\s*\|/.test(line) && /^\s*\|[\s:|-]+\|\s*$/.test(lines[i + 1] || '')) {
      const rows = [];
      while (i < lines.length && /^\s*\|/.test(lines[i])) rows.push(lines[i++]);
      out.push(makeTable(rows));
      out.push(new Paragraph({ text: '', spacing: { after: 180 } }));
      continue;
    }

    // headings
    let h = line.match(/^(#{1,4})\s+(.*)$/);
    if (h) {
      const lvl = h[1].length, txt = h[2];
      const map = { 1: HeadingLevel.HEADING_1, 2: HeadingLevel.HEADING_2, 3: HeadingLevel.HEADING_3, 4: HeadingLevel.HEADING_4 };
      const sizes = { 1: 32, 2: 27, 3: 25, 4: 24 };
      out.push(new Paragraph({
        heading: map[lvl],
        spacing: { before: lvl === 1 ? 400 : 320, after: 160, line: 300 },
        children: inline(txt, { bold: true, size: sizes[lvl], color: '000000' }),
      }));
      i++; continue;
    }

    // horizontal rule
    if (/^---+$/.test(line.trim())) {
      out.push(new Paragraph({
        text: '',
        border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: 'BBBBBB', space: 1 } },
        spacing: { before: 120, after: 200 },
      }));
      i++; continue;
    }

    // blockquote
    if (/^>\s?/.test(line)) {
      const buf = [];
      while (i < lines.length && /^>\s?/.test(lines[i])) buf.push(lines[i++].replace(/^>\s?/, ''));
      out.push(para(buf.join(' ').trim(), {
        indent: { left: 480 }, run: { italics: true }, line: 300, after: 200,
      }));
      continue;
    }

    // bullets
    if (/^[-*]\s+/.test(line)) {
      while (i < lines.length && /^[-*]\s+/.test(lines[i])) {
        out.push(new Paragraph({
          children: inline(lines[i].replace(/^[-*]\s+/, '')),
          bullet: { level: 0 },
          spacing: { line: 320, after: 80 },
        }));
        i++;
      }
      continue;
    }

    // ordered list
    if (/^\d+\.\s+/.test(line)) {
      while (i < lines.length && /^\d+\.\s+/.test(lines[i])) {
        out.push(new Paragraph({
          children: inline(lines[i]),
          indent: { left: 360, hanging: 180 },
          spacing: { line: 320, after: 80 },
        }));
        i++;
      }
      continue;
    }

    // paragraph
    const buf = [];
    while (i < lines.length && lines[i].trim() &&
           !/^(#{1,4})\s/.test(lines[i]) && !/^>/.test(lines[i]) &&
           !/^\s*\|/.test(lines[i]) && !/^[-*]\s+/.test(lines[i]) &&
           !/^\d+\.\s+/.test(lines[i]) && !/^---+$/.test(lines[i].trim())) {
      buf.push(lines[i++]);
    }
    out.push(para(buf.join(' ').trim()));
  }
  return out;
}

// ---- assemble ---------------------------------------------------------------
const files = [
  'manuscript/00-abstract.md',
  'manuscript/01-introduction.md',
  'manuscript/02-context.md',
  'manuscript/03-approach.md',
  'manuscript/04-findings.md',
  'manuscript/05-conclusion.md',
  'manuscript/06-references.md',
];

const children = [];
files.forEach((f, idx) => {
  let md = fs.readFileSync(f, 'utf8');

  // strip the per-section reference appendix from the methodology
  // split the author notes off the references file and mark them
  if (f.includes('06-references')) {
    const parts = md.split('## Notes for the author before submission');
    md = parts[0].replace(/---\s*$/, '');
    if (parts[1]) {
      children.push(...parse(md));
      children.push(new Paragraph({ children: [new PageBreak()] }));
      children.push(new Paragraph({
        children: inline('APPENDIX — EDITORIAL NOTES (DELETE BEFORE SUBMISSION)', { bold: true, size: 28, color: '990000' }),
        spacing: { before: 200, after: 200 },
      }));
      children.push(...parse('## Notes for the author before submission' + parts[1]));
      return;
    }
  }

  if (idx > 0) children.push(new Paragraph({ children: [new PageBreak()] }));
  children.push(...parse(md));
});

const doc = new Document({
  creator: 'Samet Eker',
  title: "Whose 'we'? How Türkiye's Maarif Model and the International Baccalaureate position the student",
  styles: {
    default: {
      document: { run: { font: FONT, size: SIZE }, paragraph: { spacing: { line: 360 } } },
    },
  },
  sections: [{
    properties: { page: { margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    children,
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync('Whose-we-Maarif-Model-and-IB-learner-profile.docx', buf);
  console.log('written:', buf.length, 'bytes');
});
