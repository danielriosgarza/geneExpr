#!/usr/bin/env python3
"""
Generate a single HTML report with interactive Plotly volcano plots for each
comparison group across three species. WebGL traces are converted to SVG so
markers remain visible.
"""
import os
import numpy as np
import cobra
import plotly.graph_objects as go
from parseGenExpData import GeneExpr, analyze_group


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_PATHS = {
    'bh': '/users/drgarza/documents/github/geneexpr/files/bh',
    'bt': '/users/drgarza/documents/github/geneexpr/files/bt',
    'ri': '/users/drgarza/documents/github/geneexpr/files/ri',
}

DISPLAY_NAMES = {
    'bh': 'Blautia hydrogenotrophica',
    'bt': 'Bacteroides thetaiotaomicron',
    'ri': 'Roseburia intestinalis',
}

GROUP_SETTINGS = {
    'bh': {
        'rootID': 'ncbi',
        'labels': ['t14', 't32', 't72'],
        'idxs': [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
        'comparisons': {
            ('14h', '32h'): 't14vst32_deseq.txt',
            ('14h', '72h'): 't14vst72_deseq.txt',
            ('32h', '72h'): 't32vst72_deseq.txt',
        },
    },
    'bt': {
        'rootID': None,
        'labels': ['t04', 't12', 't36'],
        'idxs': [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
        'comparisons': {
            ('04h', '12h'): 't04vst12_deseq.txt',
            ('04h', '36h'): 't04vst36_deseq.txt',
            ('12h', '36h'): 't12vst36_deseq.txt',
        },
    },
    'ri': {
        'rootID': None,
        'labels': ['t04', 't12', 't48'],
        'idxs': [[4, 5, 6], [7, 8, 9], [10, 11, 12]],
        'comparisons': {
            ('04h', '12h'): 't04vst12_deseq.txt',
            ('04h', '48h'): 't04vst48_deseq.txt',
            ('12h', '48h'): 't12vst48_deseq.txt',
        },
    },
}

# -----------------------------------------------------------------------------
# Helper to convert scattergl → scatter (SVG)
# -----------------------------------------------------------------------------

def to_svg(fig: go.Figure) -> go.Figure:
    svg = go.Figure(layout=fig.layout)
    for tr in fig.data:
        if tr.type == 'scattergl':
            marker_dict = (
                tr.marker.to_plotly_json() if hasattr(tr.marker, 'to_plotly_json')
                else dict(tr.marker)
            )
            svg.add_trace(
                go.Scatter(
                    x=tr.x,
                    y=tr.y,
                    mode=tr.mode,
                    marker=marker_dict,
                    name=tr.name,
                    customdata=getattr(tr, 'customdata', None),
                    hovertemplate=getattr(tr, 'hovertemplate', None),
                )
            )
        else:
            svg.add_trace(tr)
    return svg


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    html_sections = []

    for sp, base_dir in BASE_PATHS.items():
        gene_dir = os.path.join(base_dir, 'genes')
        model_path = os.path.join(base_dir, 'model', f'{sp}_final.xml')
        wc_path = os.path.join(gene_dir, 'wcReactions.txt')

        # Load model
        model = cobra.io.read_sbml_model(model_path)
        settings = GROUP_SETTINGS[sp]

        # Instantiate GeneExpr
        ge = GeneExpr(
            geneFolder=gene_dir,
            tpmFile=f'{sp}_tpm.txt',
            groupLabels=settings['labels'],
            groupIDX=settings['idxs'],
            groupComparison=settings['comparisons'],
            featuresFile=f'{sp}_BVBRC_features.txt',
            sbmlModel=model,
            rootID=settings['rootID'],
            species=sp,
        )

        reaction_ids = np.array([r.id for r in model.reactions])
        with open(wc_path) as f:
            wc_reactions = [ln.strip() for ln in f if ln.strip()]

        # Build volcano for each comparison
        for gA, gB in settings['comparisons']:
            fig, _ = analyze_group(
                gA, gB, ge, model, reaction_ids, wc_reactions, DISPLAY_NAMES[sp]
            )
            svg_fig = to_svg(fig)
            div = svg_fig.to_html(full_html=False, include_plotlyjs='cdn')
            html_sections.append(
                f"<h2>{DISPLAY_NAMES[sp]}: {gA} vs {gB}</h2>\n{div}\n"
            )

    # Write report
    report_html = (
        '<html><head>'
        '<meta charset="utf-8">'
        '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
        '<title>Gene Expression Comparisons</title>'
        '</head><body>\n' + '\n'.join(html_sections) + '</body></html>'
    )

    with open('comparison_report.html', 'w') as fh:
        fh.write(report_html)
    print('Generated interactive report: comparison_report.html')


if __name__ == '__main__':
    main()