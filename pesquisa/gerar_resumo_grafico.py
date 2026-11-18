#!/usr/bin/env python3
"""
Resumo gráfico (graphical abstract) do artigo (17/11/2026, diário Entrada 153).
Exigência do periódico: imagem única, mínimo de 531 x 1328 px (altura x largura), legível em 5 x 13 cm; TIFF, EPS, PDF ou Office.
Desenhado em 13 x 5,2 cm (as fontes de 6 a 15 pt valem nesse tamanho) e gravado em PDF (vetorial) e PNG a 600 dpi (3071 x 1228 px).
Sem identificação dos autores. Cores neutras (a rampa de cinza-azulado das Figuras S1 e S2) e um único acento para os números.
Todos os números vêm do próprio artigo (resumo e Seção 6): 0,011 p.p. (regret máximo, Tabela S19), 44% (Tabela 2), 99-100% (Tabela S18),
6-9% (Tabela 7, 5,7% e 8,8%). O script confere que nenhum texto passa da borda da sua caixa.
Saída: pesquisa/latex_ejor_c/figures/graphical_abstract.{pdf,png}
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch   # noqa: E402

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(BASE, 'pesquisa', 'latex_ejor_c', 'figures')
os.makedirs(OUT, exist_ok=True)

INK, INK2 = '#1f252d', '#4b5361'
N2, N3 = '#69727f', '#323944'                     # da rampa neutra das Figuras S1 e S2
FILL, LINE = '#eef0f3', '#c9cdd4'
ACENTO = '#2b6a8f'                                # um único acento, para o número que importa (contraste 5,7:1 sobre branco)
W, H = 130.0, 52.0                                # mm
MARG, GAP = 0.8, 6.6
PW = (W - 2 * MARG - 2 * GAP) / 3                 # largura de cada painel
X = [MARG + k * (PW + GAP) for k in range(3)]

plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'figure.dpi': 600})
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
CONFERIR = []     # (texto, x_min, x_max) em mm


def caixa(x, y, w, h, fc=FILL, ec=LINE, lw=0.5, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw))


def texto(x, y, s, lim, **kw):
    t = ax.text(x, y, s, **kw); CONFERIR.append((t, lim[0], lim[1])); return t


def numero(px, num, y_num, y_txt, legenda, fs_num=15.0, fs_txt=6.2):
    cx = px + PW / 2
    texto(cx, y_num, num, (px + 1.2, px + PW - 1.2), ha='center', va='center', fontsize=fs_num, fontweight='bold', color=ACENTO)
    texto(cx, y_txt, legenda, (px + 1.2, px + PW - 1.2), ha='center', va='center', fontsize=fs_txt, color=INK2, linespacing=1.2)


def cabecalho(px, n, s):
    ax.add_patch(plt.Circle((px + 2.3, 47.4), 2.1, fc=N3, ec='none'))
    ax.text(px + 2.3, 47.35, str(n), ha='center', va='center', fontsize=7.4, fontweight='bold', color='white')
    texto(px + 5.6, 47.4, s, (px + 5.6, px + PW), ha='left', va='center', fontsize=7.8, fontweight='bold', color=INK)


# ── 1. o que importa ────────────────────────────────────────────────────────
cabecalho(X[0], 1, 'What matters')
caixa(X[0], 24.2, PW, 20.2); caixa(X[0], 2.0, PW, 20.2)
numero(X[0], '0.011 p.p.', 38.4, 29.6, 'Quality ranking: max regret\nover 23 thermal and\ntravel-time conditions')
numero(X[0], '44%', 16.2, 7.4, 'Punctuality: high-uncertainty\nruns of deterministic plans\nwith no on-time solution')

# ── 2. o fluxo ──────────────────────────────────────────────────────────────
cabecalho(X[1], 2, 'Workflow')
passos = ['Inflate travel times to\na calibrated quantile', 'Solve with any routing\nsolver, even a strong one',
          'Certify by simulation:\non time with P ≥ 0.90', 'Refine the plan if\npunctuality still holds']
alt, esp = 8.6, 2.6
for k, s in enumerate(passos):
    y = 44.4 - alt - k * (alt + esp)
    caixa(X[1], y, PW, alt, fc='#ffffff', ec=N2, lw=0.7)
    ax.add_patch(plt.Circle((X[1] + 3.4, y + alt / 2), 1.7, fc=N2, ec='none'))
    ax.text(X[1] + 3.4, y + alt / 2 - 0.05, chr(ord('a') + k), ha='center', va='center', fontsize=6.6, fontweight='bold', color='white')
    texto(X[1] + 6.7, y + alt / 2, s, (X[1] + 6.7, X[1] + PW - 0.8), ha='left', va='center', fontsize=6.4, color=INK, linespacing=1.15)
    if k < len(passos) - 1:
        ax.add_patch(FancyArrowPatch((X[1] + PW / 2, y - 0.15), (X[1] + PW / 2, y - esp + 0.15), arrowstyle='-|>', mutation_scale=4.6,
                                     lw=0.8, color=N2, shrinkA=0, shrinkB=0))

# ── 3. o que se ganha ───────────────────────────────────────────────────────
cabecalho(X[2], 3, 'What it buys')
caixa(X[2], 24.2, PW, 20.2); caixa(X[2], 2.0, PW, 20.2)
numero(X[2], '99–100%', 38.4, 29.6, 'of runs comply out of\nsample (lower-bound rule)')
numero(X[2], '6–9%', 16.2, 7.4, 'lower modeled cost with a\nstrong solver, 200–400\ncustomers (mostly the solver)')

# setas entre os painéis
for k in (0, 1):
    x0 = X[k] + PW + 0.9
    ax.add_patch(FancyArrowPatch((x0, 24.0), (x0 + GAP - 1.8, 24.0), arrowstyle='-|>', mutation_scale=7, lw=1.1, color=N2, shrinkA=0, shrinkB=0))

# confere que cada texto cabe na sua caixa (em mm)
fig.canvas.draw()
rend = fig.canvas.get_renderer()
inv = ax.transData.inverted()
erros = 0
for t, lo, hi in CONFERIR:
    bb = t.get_window_extent(rend)
    (x0, _), (x1, _) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1))
    if x0 < lo - 0.05 or x1 > hi + 0.05:
        erros += 1
        print(f'  TEXTO FORA DA CAIXA: {t.get_text()!r}: {x0:.1f}-{x1:.1f} mm, limite {lo:.1f}-{hi:.1f} mm')
print('conferência de texto nas caixas:', 'OK' if not erros else f'{erros} texto(s) passam da caixa')

fig.savefig(os.path.join(OUT, 'graphical_abstract.pdf'), facecolor='white')
fig.savefig(os.path.join(OUT, 'graphical_abstract.png'), dpi=600, facecolor='white')
print('resumo gráfico:', sorted(f for f in os.listdir(OUT) if f.startswith('graphical')))
