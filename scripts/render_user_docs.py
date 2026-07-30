#!/usr/bin/env python3
"""Render Markdown user guides under docs/user/ to HTML next to pdoc output.

Used by ``make docu`` so installation and other how-tos ship with the API site.
Stdlib only (no extra Markdown dependency): supports a practical subset of GFM.
"""

from __future__ import annotations

import html
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
USER_MD_DIR = REPO_ROOT / "docs" / "user"
DOCS_OUT = REPO_ROOT / "docs"

# Source markdown → output HTML basename
PAGES = {
    "installation.md": "installation.html",
}


def _slugify(text: str) -> str:
    s = text.strip().lower()
    s = re.sub(r"[^\w\s-]", "", s)
    s = re.sub(r"[-\s]+", "-", s)
    return s.strip("-")


def md_to_html_body(md: str) -> str:
    """Convert a constrained Markdown subset to HTML body fragments."""
    lines = md.replace("\r\n", "\n").split("\n")
    out: list[str] = []
    i = 0
    in_code = False
    code_lang = ""
    code_buf: list[str] = []
    in_ul = False
    in_table = False
    table_header: list[str] | None = None

    def close_ul():
        nonlocal in_ul
        if in_ul:
            out.append("</ul>")
            in_ul = False

    def close_table():
        nonlocal in_table, table_header
        if in_table:
            out.append("</tbody></table>")
            in_table = False
            table_header = None

    def flush_para(buf: list[str]):
        if not buf:
            return
        text = " ".join(buf).strip()
        if text:
            out.append(f"<p>{inline(text)}</p>")
        buf.clear()

    def inline(text: str) -> str:
        # Escape first, then apply simple patterns on escaped text carefully
        # Work on raw then escape segments — use ordered transforms on raw
        parts: list[str] = []
        pos = 0
        # code `...`
        pattern = re.compile(
            r"(`[^`]+`|\*\*[^*]+\*\*|\[[^\]]+\]\([^)]+\))"
        )
        for m in pattern.finditer(text):
            if m.start() > pos:
                parts.append(html.escape(text[pos : m.start()]))
            token = m.group(0)
            if token.startswith("`"):
                parts.append(f"<code>{html.escape(token[1:-1])}</code>")
            elif token.startswith("**"):
                parts.append(f"<strong>{html.escape(token[2:-2])}</strong>")
            elif token.startswith("["):
                label, url = re.match(r"\[([^\]]+)\]\(([^)]+)\)", token).groups()
                # rewrite relative md links to html counterparts where known
                href = url
                if url.endswith("installation.md") or url.endswith("./installation.md"):
                    href = "installation.html"
                elif "Installing" in url:
                    href = "installation.html"
                parts.append(
                    f'<a href="{html.escape(href, quote=True)}">{html.escape(label)}</a>'
                )
            pos = m.end()
        if pos < len(text):
            parts.append(html.escape(text[pos:]))
        return "".join(parts)

    para: list[str] = []

    while i < len(lines):
        line = lines[i]

        if line.startswith("```"):
            close_ul()
            close_table()
            flush_para(para)
            if not in_code:
                in_code = True
                code_lang = line[3:].strip()
                code_buf = []
            else:
                lang_attr = f' class="language-{html.escape(code_lang)}"' if code_lang else ""
                code_html = html.escape("\n".join(code_buf))
                out.append(f"<pre><code{lang_attr}>{code_html}</code></pre>")
                in_code = False
            i += 1
            continue

        if in_code:
            code_buf.append(line)
            i += 1
            continue

        # table rows
        if "|" in line and line.strip().startswith("|"):
            close_ul()
            flush_para(para)
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            # separator row
            if all(re.match(r"^:?-+:?$", c or "") for c in cells):
                i += 1
                continue
            if not in_table:
                out.append('<table class="user-doc-table">')
                out.append("<thead><tr>")
                for c in cells:
                    out.append(f"<th>{inline(c)}</th>")
                out.append("</tr></thead><tbody>")
                in_table = True
                table_header = cells
            else:
                out.append("<tr>")
                for c in cells:
                    out.append(f"<td>{inline(c)}</td>")
                out.append("</tr>")
            i += 1
            continue
        else:
            close_table()

        if not line.strip():
            close_ul()
            flush_para(para)
            i += 1
            continue

        heading = re.match(r"^(#{1,6})\s+(.*)$", line)
        if heading:
            close_ul()
            flush_para(para)
            level = len(heading.group(1))
            title = heading.group(2).strip()
            slug = _slugify(title)
            out.append(
                f'<h{level} id="{html.escape(slug, quote=True)}">{inline(title)}</h{level}>'
            )
            i += 1
            continue

        if re.match(r"^[-*]\s+", line):
            flush_para(para)
            if not in_ul:
                out.append("<ul>")
                in_ul = True
            item = re.sub(r"^[-*]\s+", "", line)
            out.append(f"<li>{inline(item)}</li>")
            i += 1
            continue

        if re.match(r"^\d+\.\s+", line):
            close_ul()
            flush_para(para)
            items = [re.sub(r"^\d+\.\s+", "", line)]
            i += 1
            while i < len(lines) and re.match(r"^\d+\.\s+", lines[i]):
                items.append(re.sub(r"^\d+\.\s+", "", lines[i]))
                i += 1
            out.append("<ol>")
            for it in items:
                out.append(f"<li>{inline(it)}</li>")
            out.append("</ol>")
            continue

        close_ul()
        para.append(line.strip())
        i += 1

    close_ul()
    close_table()
    flush_para(para)
    if in_code:
        raise SystemExit("Unclosed fenced code block in Markdown")

    return "\n".join(out)


def page_shell(title: str, body: str, *, active: str = "installation") -> str:
    """Minimal HTML shell with nav links to API docs and install page."""
    nav_install = "class=\"active\"" if active == "installation" else ""
    nav_api = "class=\"active\"" if active == "api" else ""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>{html.escape(title)} — CoreMS</title>
  <style>
    :root {{
      --bg: #fafafa;
      --fg: #1a1a1a;
      --muted: #555;
      --link: #0b57d0;
      --border: #e0e0e0;
      --code-bg: #f0f0f0;
      --nav-bg: #111827;
      --nav-fg: #f9fafb;
    }}
    @media (prefers-color-scheme: dark) {{
      :root {{
        --bg: #0f1115;
        --fg: #e8eaed;
        --muted: #9aa0a6;
        --link: #8ab4f8;
        --border: #2a2f3a;
        --code-bg: #1c212b;
        --nav-bg: #0b0d12;
        --nav-fg: #e8eaed;
      }}
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, sans-serif;
      line-height: 1.55;
      color: var(--fg);
      background: var(--bg);
    }}
    header.site-nav {{
      background: var(--nav-bg);
      color: var(--nav-fg);
      padding: 0.75rem 1.25rem;
      display: flex;
      gap: 1.25rem;
      align-items: center;
      flex-wrap: wrap;
    }}
    header.site-nav a {{
      color: var(--nav-fg);
      text-decoration: none;
      opacity: 0.9;
    }}
    header.site-nav a:hover {{ opacity: 1; text-decoration: underline; }}
    header.site-nav a.active {{ font-weight: 700; text-decoration: underline; }}
    header.site-nav .brand {{ font-weight: 700; margin-right: 0.5rem; }}
    main {{
      max-width: 52rem;
      margin: 0 auto;
      padding: 1.5rem 1.25rem 3rem;
    }}
    h1, h2, h3, h4 {{ line-height: 1.25; }}
    h1 {{ font-size: 1.85rem; margin-top: 0; }}
    h2 {{ font-size: 1.35rem; margin-top: 2rem; border-bottom: 1px solid var(--border); padding-bottom: 0.25rem; }}
    h3 {{ font-size: 1.1rem; margin-top: 1.4rem; }}
    a {{ color: var(--link); }}
    code {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 0.92em;
      background: var(--code-bg);
      padding: 0.1em 0.35em;
      border-radius: 4px;
    }}
    pre {{
      background: var(--code-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 0.9rem 1rem;
      overflow-x: auto;
    }}
    pre code {{ background: transparent; padding: 0; }}
    table.user-doc-table {{
      border-collapse: collapse;
      width: 100%;
      margin: 1rem 0;
      font-size: 0.95rem;
    }}
    table.user-doc-table th, table.user-doc-table td {{
      border: 1px solid var(--border);
      padding: 0.45rem 0.6rem;
      text-align: left;
      vertical-align: top;
    }}
    table.user-doc-table th {{ background: var(--code-bg); }}
    ul, ol {{ padding-left: 1.4rem; }}
    p {{ margin: 0.75rem 0; }}
    footer {{
      max-width: 52rem;
      margin: 0 auto;
      padding: 0 1.25rem 2rem;
      color: var(--muted);
      font-size: 0.9rem;
    }}
  </style>
</head>
<body>
  <header class="site-nav">
    <span class="brand">CoreMS</span>
    <a href="installation.html" {nav_install}>Installation</a>
    <a href="corems.html" {nav_api}>API reference</a>
    <a href="https://github.com/EMSL-Computing/CoreMS">GitHub</a>
  </header>
  <main>
{body}
  </main>
  <footer>
    Generated from <code>docs/user/</code> Markdown via <code>make docu</code> / <code>scripts/render_user_docs.py</code>.
  </footer>
</body>
</html>
"""


def main() -> None:
    if not USER_MD_DIR.is_dir():
        raise SystemExit(f"Missing user docs directory: {USER_MD_DIR}")

    DOCS_OUT.mkdir(parents=True, exist_ok=True)

    for src_name, out_name in PAGES.items():
        src = USER_MD_DIR / src_name
        if not src.is_file():
            raise SystemExit(f"Missing user doc: {src}")
        md = src.read_text(encoding="utf-8")
        # Drop leading H1 for title extraction
        title = "CoreMS"
        m = re.match(r"^#\s+(.*)$", md, re.M)
        if m:
            title = m.group(1).strip()
        body = md_to_html_body(md)
        out_path = DOCS_OUT / out_name
        out_path.write_text(
            page_shell(title, body, active="installation"),
            encoding="utf-8",
        )
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}")

    # Landing index: prefer installation as default entry (API still one click away)
    index = DOCS_OUT / "index.html"
    index.write_text(
        """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>CoreMS documentation</title>
  <meta http-equiv="refresh" content="0; url=./installation.html"/>
  <link rel="canonical" href="./installation.html"/>
</head>
<body>
  <p>Redirecting to <a href="./installation.html">Installation</a>.
  API reference: <a href="./corems.html">corems.html</a>.</p>
</body>
</html>
""",
        encoding="utf-8",
    )
    print(f"Wrote {index.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
