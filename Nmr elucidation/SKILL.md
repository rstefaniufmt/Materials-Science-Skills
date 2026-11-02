---
name: nmr-elucidation
description: "Lê, processa e interpreta RMN 1D de ¹H e ¹³C (com DEPT-135/90 e APT) a partir de JCAMP-DX (.jdx/.dx, espectro ou FID; Bruker, JEOL, Varian, MestReNova, Nanalysis, Magritek, Oxford), pastas Bruker/Varian e CSV: referência pelo solvente, fase e linha de base, picos, integrais em nº de H, multiplicidade e constantes J (árvore de acoplamentos), solventes/impurezas comuns, listas no formato de artigo, multiplicidade dos carbonos por DEPT, e ranqueamento de estruturas candidatas por predição de deslocamentos (CSP5) com atribuição átomo a átomo e probabilidade relativa estilo DP4 (desreplicação, isômeros, regioquímica). Use sempre que o usuário trouxer espectros de RMN ou pedir para atribuir sinais, achar J, montar a lista de dados ¹H/¹³C, conferir pureza/solvente residual, comparar com a literatura, decidir entre estruturas ou propor caminhos de elucidação - mesmo sem dizer 'RMN' (ex.: 'NMR do composto', 'espectro de próton', 'carbono 13', 'DEPT', 'qual é a estrutura')."
---

# RMN ¹H/¹³C: processamento, atribuição, J e estruturas candidatas

Scripts em `scripts/` (numpy, scipy, pandas, matplotlib, **nmrglue**; **rdkit** e **csp5** para
estruturas). Instale o que faltar: `pip install nmrglue rdkit csp5` (o CSP5 traz PyTorch e os modelos;
roda offline). Rode de qualquer pasta:

```bash
S=<caminho-desta-skill>/scripts
python $S/nmr.py info *.jdx                                    # núcleo, campo, solvente, referência, fase, S/N
python $S/nmr.py h1 1h.jdx --formula C9H8O4 -o h1.png --multiplets mult.png --csv h1.csv
python $S/nmr.py c13 13c.jdx --dept135 dept135.jdx --dept90 dept90.jdx -o c13.png
python $S/nmr.py report --h1 1h.jdx --c13 13c.jdx --dept135 d135.jdx --dept90 d90.jdx --formula C11H14N2O -d rel
python $S/nmr.py candidates --h1 1h.jdx --c13 13c.jdx --dept135 d135.jdx --formula C11H14N2O \
       --smiles "citisina=O=C1C=CC=C2N1CC3CNCC2C3" "isomero=..." -d cand      # ou --file candidatas.smi
python $S/nmr.py predict "CC(=O)Oc1ccccc1C(=O)O" --solvent cdcl3 -o previsto.png
python $S/nmr.py impurities --solvent dmso-d6                  # Gottlieb/Fulmer
python $S/nmr.py process fid.jdx -o espectro.jdx                # FID -> espectro (fase automática)
```

Opções comuns: `--solvent` (se o arquivo não diz), `--ref auto|file|tms|7.31=7.26`, `--norm 2.29=3`
(fixa uma integral), `--formula` (normaliza integrais pelo total de H e checa candidatas), `-k` (limiar
de picos ¹H em σ), `--decimal ,`, `--lang en`.

Em Python: `sys.path.insert(0, S)` e use `nmr_io.load`, `nmr_process.prepare` (linha de base,
referência, refase), `nmr_peaks` (pick_peaks, analyze_1h, analyze_13c, acs_1h, acs_13c, inventory),
`nmr_struct` (Candidate, rank, mol_info, predict, draw_assigned), `nmr_data` (solventes, impurezas,
faixas de δ, J típicos, IDH) e `nmr_plot`.

## O que esta skill faz - e o que não faz

Faz bem: ler/processar, **atribuir sinais e extrair J de multipletos de 1ª ordem**, contar H e C
(com DEPT), apontar solvente/água/impurezas/satélites, montar o texto da parte experimental, resumir
o que os dados dizem (IDH, H aromáticos, CH₃ isolados, etilas...), e **escolher entre estruturas
candidatas** (desreplicação de compostos conhecidos, regioisômeros, glicosídeos, confirmação de
síntese) com atribuição átomo a átomo.

Não faz: **elucidar do zero uma estrutura inédita** só com ¹H/¹³C 1D. Isso exige fórmula (EM de alta
resolução), conectividade (HSQC, HMBC, COSY; NOESY/ROESY para estereoquímica) e geração de estruturas
(CASE) - e mesmo assim termina em candidatas que se ranqueiam como aqui. A probabilidade relativa só vale
**dentro** da lista dada: se a estrutura certa não estiver nela, a "vencedora" é só a menos errada.
Diga isso ao usuário sempre que ranquear.

## Fluxo de trabalho

1. **`info` primeiro.** Confira núcleo, campo, solvente e, principalmente, a **referência**: a skill
   reconhece o sinal residual do solvente pelo padrão de acoplamento com ²H (CDCl3 ¹³C tripleto 1:1:1,
   DMSO-d6 septeto, CHD₂ quintetos) e corrige desvios pequenos (≤ 0,1 ppm ¹H, ≤ 1 ppm ¹³C); desvio maior
   é só avisado (`--ref obs=esperado`). Fase ruim em espectro processado é refeita (Hilbert + ACME). FID é
   transformado (filtro digital Bruker, LB 0,3/1 Hz, zero-filling, fase automática).
2. **¹H**: `h1` com `--formula` sempre que houver (normaliza integrais pelo total de H e avisa quando
   faltam H trocáveis). Leia a tabela: `mult`/`J` (dd, td, ddd... por desconvolução de 1ª ordem; `m` quando
   não há solução), `nota` (telhado/2ª ordem, sobreposição com solvente, sinal largo/trocável),
   `parceiros_J` (multipletos que compartilham J → candidatos a acoplados entre si - confirme com COSY).
   Olhe `--multiplets` (zoom com árvore de J) antes de reportar J.
3. **¹³C**: `c13`, com DEPT-135/90 (ou APT) se houver → CH₃/CH₂/CH/C e nº de H em carbono (a diferença
   para o total de H da fórmula = H em heteroátomos). Picos minúsculos, solvente e impurezas saem da lista.
4. **Interpretar**: `report` junta tudo e dá o "Resumo estrutural". Use `references/interpretation.md`
   (roteiro de elucidação, o que cada dado restringe, quando pedir 2D) e `references/shifts_couplings.md`.
5. **Candidatas**: reúna estruturas plausíveis (síntese esperada e subprodutos, isômeros de posição,
   compostos já descritos na espécie/gênero, PubChem/COCONUT/LOTUS/NPAtlas pela fórmula) e rode
   `candidates`. Leia `references/candidates.md` antes de concluir: MAE ¹³C ≲ 1,5-2 ppm e ¹H ≲ 0,15-0,2 ppm
   com DEPT e integrais coerentes = compatível; prob. relativa alta com MAE ruim = nenhuma serve.
6. **Entregue**: lista de dados no formato de artigo, tabela de atribuição (δC, δH, mult., J) e, se
   ranqueou, o veredito com números e a ressalva de que é relativo à lista. Mostre as figuras (abra antes).

## Regras que evitam erros comuns

- **Integrais de ¹H** dependem de d1 ≥ 5·T1 (aromáticos/CH₃ isolados relaxam devagar) e de linha de base;
  diferenças de ±0,1-0,2 H são normais. **¹³C desacoplado não é quantitativo** (NOE, T1 de C quaternário).
- **H trocáveis** (OH, NH, COOH) somem em CD3OD/D₂O, alargam e mudam de lugar com concentração/água:
  não force atribuição; confirme por troca com D₂O.
- **J só de multipletos de 1ª ordem** (Δν/J ≳ 5). Telhado forte, campo baixo (60-90 MHz) e sistemas
  AA'BB' dão "dubletos" falsos: a skill avisa com `nota`; reporte como m quando duvidoso.
- **Mesmo J em dois sinais** sugere acoplamento entre eles, mas J típicos (7-8 Hz) se repetem por acaso.
- **Sinal no lugar de impureza** (acetona 2,17; água 1,56; graxa 1,26/0,86 em CDCl3) só é impureza se
  for pequeno e/ou variar entre amostras; a skill só rotula sinais pequenos.
- **Simetria**: menos sinais de ¹³C que carbonos = simetria (ou sobreposição); mais = mistura, rotâmeros
  (amidas, carbamatos), diastereoisômeros ou impureza.
- **Formato**: δ com 2 casas (¹H) e 1 casa (¹³C), J com 1 casa, multiplicidades s, d, t, q, quint, sext,
  sept, dd, dt, td, ddd, m, sl (sinal largo / br s).

## Formatos

| Origem | Lido | Observações |
|---|---|---|
| JCAMP-DX 4.24/5/6 (.jdx, .dx) espectro (XYDATA ou NTUPLES) | sim | Bruker TopSpin, MestReNova, JEOL Delta, VnmrJ, NMRium, Spinsolve, Nanalysis; eixo em Hz ou ppm |
| JCAMP-DX FID (NMR FID) | sim | filtro digital Bruker removido; LB/zero-fill/fase automáticos (`--lb`, `process`) |
| Pasta Bruker (fid + acqus ou pdata/1/1r) e Varian/Agilent (fid + procpar) | sim (nmrglue) | |
| CSV/TXT (ppm, intensidade) ou lista de picos | sim | informe `--sf` (MHz) para J em Hz |
| .mnova, .jdf (JEOL), .1r isolado sem procs | não | exporte JCAMP-DX (ver `references/formats.md`) |

## Arquivos de referência

- `references/interpretation.md` - roteiro de elucidação com 1D, o que cada dado restringe, quando e
  quais 2D pedir, desreplicação, como apresentar a atribuição.
- `references/shifts_couplings.md` - faixas de δ ¹H/¹³C, J típicos, efeitos de substituintes, sinais
  diagnósticos (açúcares, flavonoides, terpenos, ésteres, amidas).
- `references/candidates.md` - predição CSP5, atribuição ótima, escores, prob. estilo DP4 e seus limites.
- `references/processing.md` - referência, fase, linha de base, integração, DEPT/APT, satélites,
  solventes e impurezas, armadilhas de aquisição.
- `references/formats.md` - JCAMP-DX (dialetos), como exportar de TopSpin, MestReNova, Delta, VnmrJ,
  Spinsolve; depuração da leitura.
