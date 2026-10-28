# PROMPT PARA NOVA CONVERSA — VirtualTissue-scRNAseq / melanoma + anti-PD-1

> Copie o conteúdo deste arquivo para uma nova conversa, ou peça ao agente que
> leia este arquivo inteiro e execute suas instruções. O texto não depende do
> histórico do chat. Inventário local conferido em 28/10/2026.

## 1. Meu pedido: implemente o produto e entregue uma versão testável

Sou Helder Nakaya. Quero que você implemente e conclua uma versão de
**VirtualTissue-scRNAseq para simular melanoma sob anti-PD-1 a partir de células
obtidas por 10x single-cell RNA-seq de tumor dissociado**, pronta para eu subir no
GitHub e testar com meu grupo.

Quero carregar uma biópsia, gerar um tecido virtual e acompanhar células
individuais: migração, encontros, reconhecimento antigênico assumido ou
informado, ativação, diferenciação/transições de estado, disfunção/exaustão,
proliferação, citotoxicidade, sobrevivência e morte. Quero inspecionar por que
uma célula permaneceu em um estado ou mudou, com probabilidades/hazards,
condições locais, exposição ao anti-PD-1 e histórico de eventos.

**O resultado não pode ser apenas uma animação, uma grade de pontos ou células
que só alternam entre mover/atacar/descansar.** As transições de estado de CD8
precisam existir no motor e alterar capacidades subsequentes. Fadiga numérica
após atacar não é sinônimo de exaustão. Estados de CD8 não são tipos celulares
amplos. Não converta CD8 em outra linhagem ao chamar algo de diferenciação.

O software deve funcionar com os dados próprios do usuário, com um exemplo real
público de melanoma e com um pequeno exemplo sintético offline para testes.
Use como referência de experiência, inspeção e qualidade o projeto
https://github.com/hylixinsights/VirtualTissue-Xenium .

Execute o trabalho; não termine apenas com uma proposta, diagnóstico ou outro
prompt. Reaproveite o código útil, corrija o que está errado e persista até
entregar os critérios de aceitação abaixo. Comunique-se comigo em português;
código, documentação do repositório e notebooks devem ficar em inglês.

## 2. Escopo imediato e compromisso com o plano científico

Existem dois objetivos que não devem ser confundidos:

1. **Entrega de engenharia agora:** simulador mecanístico exploratório, com
   estados dinâmicos de CD8, anti-PD-1 explícito, hipóteses versionadas, controles,
   dados próprios, visualizador e testes. Pode usar parâmetros assumidos quando
   não houver estimativas defensáveis; essa condição precisa aparecer no pacote,
   na interface e nos resultados. Este é o produto que quero testar no grupo.
2. **Validação do plano SPECIAL:** estimar parâmetros de biópsias pareadas e
   demonstrar previsão da composição dos estados de CD8 em pacientes novos.
   Não declarar isso alcançado sem executar os respectivos testes/gates.

A insuficiência da coorte para o desenho confirmatório **não deve paralisar a
implementação e verificação do produto exploratório**. Ela impede alegar
validação clínica/preditiva e ajustar resultados para aparentar sucesso.
Por outro lado, construir uma demonstração não completa as tarefas científicas
originais. Mantenha dois quadros de status: implementação do produto e estudo
preditivo. Quando uma tarefa científica não puder ser executada, registre o
motivo concreto e continue as tarefas de engenharia independentes.

Registre essa prioridade atual em uma **nova emenda prospectiva**, sem alterar
`prereg/FREEZE.json`, os planos congelados ou `config/study.yaml`. A emenda deve
separar implementação de modelos/infraestrutura e avaliação exploratória de
uma aprovação de G1/G4/G5. Não inspecione mudanças pré/on de SPECIAL nem use
metadados de desfecho para escolher regras da demonstração. Planeje qualquer
análise exploratória nova antes de executá-la e preserve sua rastreabilidade.

**Não envie e-mail aos autores e não espere respostas deles.** Use depósitos,
artigos, suplementos e metadados públicos. Campo não recuperável fica desconhecido;
biblioteca ausente não vira zero biológico. Não use acesso controlado como
pré-requisito quando dados processados públicos bastam. Preserve as licenças.
Não faça chamadas pagas a Jev/LLM: interface de provedor e mocks podem ser
implementados; chamadas reais exigem minha autorização explícita e teto.

Prepare commits locais e uma versão pronta para publicação. **Não faça push,
release, criação de PR ou divulgação externa automaticamente**: neste pedido
quero o material pronto para eu subir e testar. Não crie outra conversa por conta
própria. Não altere o engine validado de VirtualTissue-Xenium.

## 3. Onde trabalhar e o que ler primeiro

Caminhos absolutos nesta máquina:

```text
WORKSPACE = /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue
SC_FOLDER = /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue/scRNAseq
REPO      = /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue/scRNAseq/VirtualTissue-SC
DATA      = /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue/scRNAseq/data
XENIUM    = /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue/VirtualTissue-Xenium
```

Todos os caminhos relativos abaixo são relativos a **REPO**, salvo quando
começarem com `DATA/`. `../data` a partir de REPO é DATA. Não confunda com
`REPO/data`, que contém manifests e não as matrizes grandes.

GitHub de destino: https://github.com/hylixinsights/VirtualTissue-scRNAseq.git .
Na preparação deste prompt não havia remoto configurado no checkout local;
confira novamente. Não clone outro projeto por cima do checkout existente.
`55e4469` é o commit da auditoria de alinhamento; `ade322f` contém a demonstração
v0.3 e `902155a` a demonstração anterior. O prompt pode estar em commit posterior.
Leia `git status`, `git log` e instruções AGENTS aplicáveis antes de editar.

Leitura inicial obrigatória:

1. `docs/NEW_CONVERSATION_MELANOMA_ANTI_PD1.md` — este contrato de entrega.
2. `reports/plan-alignment-audit-2026-10-28.md` — o que falta no produto atual.
3. `reports/STATUS.md`, `reports/task-board.json`, `reports/cohort_expansion_review.md`.
4. `prereg/plan-v1-original-2026-10-27.md`,
   `prereg/plan-v1.1-frozen-2026-10-27.md`,
   `prereg/working-plan-v1.2-2026-10-27.md` e `prereg/FREEZE.json`.
5. `config/study.yaml`, `config/cohort-expansion-audit.yaml`, `DATA_LICENSE.md`.
6. `data/MANIFEST.json`, `data/COHORT_EXPANSION_MANIFEST.json`,
   `reports/handoff-inventory-2026-10-28.json` e `reports/EVIDENCE_SHA256.json`.
7. `README.md`, `docs/03_model_and_assumptions.md`, `src/vtsc/tissue.py`,
   `src/vtsc/tissue_inputs.py`, `src/vtsc/model.py` e testes.
8. README e arquivos relevantes em XENIUM como referência, sem modificá-los.

A referência antiga `claude/virtualtissue-sc-special-plan-v1-2026-10-27.md`
não é o caminho de retomada deste checkout. As cópias válidas estão em `prereg/`
e em `SC_FOLDER/virtualtissue-sc-special-plan-v1.md`.

## 4. O que já existe e o que não funciona para meu objetivo

O pacote atual é v0.3.0. Há CLI `vtsc`, importação H5AD/Cell Ranger H5/MTX,
bundles com hashes, reprodução offline, notebooks e visualizador com células,
núcleos, contatos e inspetor. Há movimento, exclusão de volume, crescimento/
divisão de melanoma, dano, fadiga e remoção de cadáveres. Esses componentes
podem ser reaproveitados, com os testes necessários.

O que falta é central:

- O engine espacial mantém `cell_type` e programas iniciais fixos. Não tem
  transições Memory/Activated/Cytotoxic/Dysfunctional nem tempo-no-estado.
- A política é uma prioridade programada de ações. Não é Jev e não é uma
  política de transição de estados estimada.
- Não há exposição/schedule anti-PD-1 nem mecanismo PD-1/PD-L1 explícito.
- `recognition-off` desliga ataques; **não** representa ausência de anti-PD-1.
- `model.py` tem um processo Markoviano com taxas fornecidas, três estados
  ilustrativos e expectativa analítica. Isso não é o ajuste de M0–M3 planejado.
- O espaço não deriva do modelo não espacial. O teste de expectativa do motor
  antigo não é a equivalência ABM↔M3 exigida pelo plano.
- A anotação do exemplo 10x é uma heurística provisória, não confirmação de
  malignidade por CNV nem referência CD8 validada.

A última evidência registrada é de **34 testes e três notebooks passando**, não
validação biológica. Leia `reports/software/melanoma-validation.md` e os logs.
Não converta esses resultados em afirmações de eficácia ou previsão em pacientes.
PBMC é legado de um desvio anterior: não use como exemplo principal, nem para
substituir melanoma/TME. Não mantenha uma demonstração bonita como substituto
das transições e do tratamento pedidos aqui.

## 5. Dados reais já disponíveis, origem e processamento

### 5.1 SPECIAL: objeto imune e exportação exata

Depósito: https://doi.org/10.48804/GSAXBN , versão 8, KU Leuven RDR.
Publicações: Poźniak et al., Cell 2024 (PMID 38181739) e Roda et al., Cancer
Discovery 2025 (PMID 40530504; PMC12409281).
EGA: https://ega-archive.org/datasets/EGAD00001009291 . O raw EGA é controlado;
não é necessário pedir acesso para trabalhar com os objetos processados públicos.
Fonte pública consultável também em
https://www.biorxiv.org/content/10.1101/2023.12.14.571631.full.pdf .

Arquivo aceito:

```text
DATA/SPECIAL/raw/GC_all_immune_sharing.rds
991078732 bytes
SHA-256 93dbb75aea1ed3fd7e867c8f0f857023fdd743a316950fcbcf4a140cb35ac951
URL https://rdr.kuleuven.be/api/access/datafile/266111
```

READMEs originais: `DATA/SPECIAL/raw/Pozniak_README_v5.txt` e
`Pozniak_README_v6.txt`; URLs e hashes estão no manifest.
Licença dos dados SPECIAL/derivados: **CC BY-NC-SA 4.0**. Código original pode
continuar MIT. Não incorporar derivados SPECIAL em um pacote de dados MIT ou
produto comercial Hylix.

Exportação já concluída:

```text
DATA/SPECIAL/derived/immune/<sample_ID>/counts.mtx
DATA/SPECIAL/derived/immune/<sample_ID>/obs.parquet
DATA/SPECIAL/derived/immune/<sample_ID>/var.tsv
DATA/SPECIAL/derived/immune/<sample_ID>/author_pca.parquet
DATA/SPECIAL/derived/immune/<sample_ID>/author_harmony.parquet
DATA/SPECIAL/derived/immune/<sample_ID>/author_umap.parquet
```

São **45 pastas de amostras, 25.178 células imunes e 22.897 features**. Exemplo
de pasta existente: `sc5rCMA061`. `counts.mtx` é **genes × células**; `var.tsv`
tem cabeçalho `gene_id` e acompanha linhas; `obs.parquet` acompanha colunas e
contém `cell_id`. Não tratar isso automaticamente como diretório 10x padrão
(`matrix.mtx/features.tsv/barcodes.tsv`): implemente um adapter explícito.

O objeto original armazena versão Seurat 3.1.5. Foram preservados counts RNA,
metadados, embeddings e QC técnico recalculado; não houve reanotação/filtro.
Os valores esparsos, dimensões, somas e ordem celular passaram por round-trip
exato. Embeddings dos autores são **arquivais**, não representações independentes
para validação; metadados incluem estado/desfecho e precisam de controle de acesso.

Cadeia de scripts/evidências:

| Etapa | Script | Evidência/saída principal |
|---|---|---|
| Inventário Seurat | `r/inventory_seurat.R INPUT.rds OUTPUT.json` | `reports/evidence/immune-structure.json`; `reports/T1.1_object_inventory.md` |
| Exportação | `r/export_seurat.R INPUT.rds OUTPUT_DIR` | `DATA/SPECIAL/derived/immune/`; `reports/evidence/immune-export-validation.json` |
| Reconciliar tabelas clínicas | `scripts/extract_clinical_sources.py` | `reports/evidence/clinical-source-reconciliation.json` |
| Auditoria de pares/contagens/covariáveis | `scripts/audit_cohort.py` | `tables/samples.tsv`, `candidate_pairs.tsv`, `eligibility_pending.tsv`, `chemistry_site_timepoint.tsv`; `reports/evidence/cohort-audit.json` |
| Contagem independente no RDS | `r/validate_cohort_counts.R` | `reports/evidence/cohort-counts-r-validation.json` |
| Fluxo preliminar | `figures/plot_feasibility.py` | `figures/S1_flow_preliminary.svg` |

Não reexporte por cima das saídas aceitas. `r/export_seurat.R` exige destino
vazio e fonte registrada/verificada no manifest. Para TME/maligno, primeiro
inspecione seus schemas; adapte o leitor em vez de presumir colunas idênticas.

Tabelas clínicas locais:

- `DATA/SPECIAL/raw/clinical/Pozniak_Cell2024_TableS1.xlsx`
  — https://ars.els-cdn.com/content/image/1-s2.0-S0092867423013223-mmc1.xlsx
  — hash local registrado; checksum externo não fornecido.
- `DATA/SPECIAL/raw/clinical/cd-24-1208_supplemental_table_s1_suppst1.xlsx`
  — https://ndownloader.figshare.com/files/57656505
  — tamanho e MD5 do depósito conferidos.

Ver `data/MANIFEST.json`, `reports/T0.2_acquisition.md`,
`reports/T1.2_sample_table.md`, `reports/T1.6_confounders.md` e
`tables/S2_metadata_dictionary.tsv`.

### 5.2 SPECIAL: arquivos públicos identificados, ainda NÃO baixados

| Arquivo/destino previsto em DATA/SPECIAL/raw | URL pública | Bytes esperados | MD5 publicado |
|---|---|---:|---|
| `Entire_TME.rds` | https://rdr.kuleuven.be/api/access/datafile/210951 | 2718597537 | `fbb62e4c67f02d34c16eaa077e27d846` |
| `Malignant_cells.rds` | https://rdr.kuleuven.be/api/access/datafile/28868 | 864438613 | `f94e6421eaf685954f8c07719ebcc32a` |

Baixe e verifique esses arquivos se continuar a auditoria TME. Downloads públicos
por API funcionaram; não repetir a suposição inicial de que o navegador manual
é obrigatório. Registre tamanho, checksum e SHA-256 antes de aceitar cada arquivo.
Não baixe os ~69 GB do depósito inteiro. Os bundles Xenium SPECIAL não foram
baixados; não são requisito para executar uma geometria assumida no produto.

### 5.3 Lin2026: objeto completo baixado, somente metadados extraídos

Lin et al., Cancer Cell 2026: https://doi.org/10.1016/j.ccell.2026.08.007 .
Depósito público v2: https://doi.org/10.5281/zenodo.19586163 .

```text
DATA/Lin2026/raw/seurat_full.qs2
4021369094 bytes
MD5 a3486788163889f8cc910ea0115d70ac
SHA-256 43c9e52c5a7766fe903be29b8cbcc3aaee92c2e3a3b6801bc80a9f8b01875772
URL https://zenodo.org/records/19586163/files/seurat_full.qs2?download=1
```

Download por intervalos HTTP: `scripts/download_lin_object.py`.
`DATA/Lin2026/raw/seurat_full.qs2.part` é **incompleto, inválido para análise**;
não o confunda com o arquivo verificado. Partes em `download_parts/` são
intermediários; não redownloadar tudo sem necessidade.

Leitura integral em R/conversão oficial excedeu limites de memória. A rota bem
sucedida foi decodificar em disco e extrair apenas metadados:

```text
DATA/Lin2026/derived/audit-quarantine/source_stream.rds
  28600818752 bytes; stream intermediário, NÃO um objeto pequeno
DATA/Lin2026/derived/audit-quarantine/source_stream.integrity.json
DATA/Lin2026/derived/audit-quarantine/source_metadata.parquet
  23271965 bytes; SHA-256 e1d5b237d8166f7d8d4565bfe248dafe1740fcdce0448cf817bd9e02cb432ff6
DATA/Lin2026/derived/audit-quarantine/source_metadata.json
```

Comandos usados (não repetir sobre destinos aceitos):

```bash
.venv/bin/python scripts/decode_qs2_blocks.py ../data/Lin2026/raw/seurat_full.qs2 ../data/Lin2026/derived/audit-quarantine/source_stream.rds
.venv/bin/python scripts/extract_seurat_metadata_stream.py ../data/Lin2026/derived/audit-quarantine/source_stream.rds ../data/Lin2026/derived/audit-quarantine/source_metadata.parquet
.venv/bin/python scripts/audit_lin_object_counts.py
```

Os testes `scripts/validate_qs2_blocks.py` e `scripts/validate_metadata_stream.py`
compararam fixtures com leitores R e rejeitaram corrupção/truncamento inseguro.
Isso **não verifica a matriz real nem seu alinhamento com os metadados**. Não
usar o parquet como se já houvesse counts exportados. O objeto integrado inclui
outras coortes e reutiliza SPECIAL; não somar pacientes duplicados nem usar
embeddings integrados como referência ajustada sem vazamento.

Documentação e evidências: `docs/lin-metadata-reproduction.md`,
`reports/evidence/cohort-expansion/lin-reader-provenance.json`,
`lin-object-integrity.json`, `lin-object-count-audit.json`,
`object-count-validation.json` e `data/COHORT_EXPANSION_MANIFEST.json`.
Depósito declara CC BY 4.0; termos das coortes incorporadas continuam separados.

Fontes clínicas Lin:

- `DATA/Lin2026/raw/mmc1.pdf` —
  https://ars.els-cdn.com/content/image/1-s2.0-S1535610826003843-mmc1.pdf ;
  tabela clínica na página 11, extraída por `scripts/extract_lin_clinical.py`.
- `DATA/Lin2026/raw/mmc2.xlsx` —
  https://ars.els-cdn.com/content/image/1-s2.0-S1535610826003843-mmc2.xlsx ;
  identifica reutilização de estudos.
- `scripts/audit_cohort_expansion.py` cria os limites clínicos;
  `scripts/audit_lin_object_counts.py` cria os limites por células depositadas;
  `scripts/audit_lin_clinical_sources.py` compara fontes clínicas/EGA.
- Resultados: `tables/cohort_expansion_lin_samples.tsv`,
  `cohort_expansion_lin_object_samples.tsv`, `cohort_expansion_lin_object_patients.tsv`,
  `cohort_expansion_lin_library_crosswalk.tsv`, `cohort_expansion_reused_special.tsv`.

### 5.4 Exemplo real principal: melanoma humano dissociado da 10x

Página oficial:
https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-GEM-X .
Título: **10k Human DTC Melanoma, Chromium GEM-X Single Cell 3′**.
Licença explícita **CC BY 4.0**; um doador; sem pares pré/on ou desfecho terapêutico.
A descrição diz tecido de melanoma dissociado; o campo estruturado “Anatomical
Entity” diz “Blood”. Essa inconsistência está registrada e não deve ser ocultada.

```text
REPO/example_data/melanoma10x/filtered_feature_bc_matrix.h5
56774638 bytes; 10645 células × 38606 features
SHA-256 a2472fd15fbc708c50b95aff7f762872929a0365b9d06454cc558d29d4400091
REPO/example_data/melanoma10x/SOURCE.json
```

URL da matriz:
https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx/10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx_count_sample_filtered_feature_bc_matrix.h5

Implementação da aquisição/preparação: `src/vtsc/tissue_inputs.py`, funções
`download_melanoma`, `prepare_tissue`, `expression_features`, `annotate_demo`.
Leitor de counts: `src/vtsc/inputs.py::load_input`.
Simulação atual: `src/vtsc/tissue.py::simulate_tissue`.
Exportação/viewer: `src/vtsc/tissue_results.py`, `src/vtsc/web/tissue.html`.

Saídas existentes, todas **demonstrações de software, não tratamento anti-PD-1**:

```text
REPO/outputs/melanoma10x_input.zip
REPO/outputs/melanoma_tissue/                  # declared-recognition
REPO/outputs/melanoma_tissue.zip
REPO/outputs/melanoma_recognition_off/         # controle sem reconhecimento
REPO/outputs/melanoma_recognition_off.zip
REPO/outputs/melanoma_reopened.html
REPO/outputs/synthetic_melanoma_v03/           # fixture gerada; não paciente
```

Cada pasta de execução inclui `run.json`, `provenance.json`, `assumptions.yaml`,
`summary.json`, `decisions.jsonl`, `events.jsonl`, `trajectory.tsv`,
`manifest.json`, `explorer.html`.

O subset uniforme atual de 400 células (seed 21) contém apenas 15 provisoriamente
rotuladas melanoma-like; a fonte é rica em linfócitos. Não inventar/enriquecer
células silenciosamente para produzir uma imagem mais tumoral. Um preset
artificial de composição deve ser opcional e marcado como cenário sintético.
As contagens completas/proporcionais estão em `docs/05_example_datasets.md`.

Comandos que produziram a demonstração atual (usar nomes novos ao reproduzir):

```bash
.venv-software/bin/vtsc download-example example_data/melanoma10x
.venv-software/bin/vtsc prepare-tissue example_data/melanoma10x/filtered_feature_bc_matrix.h5 outputs/melanoma10x_input.zip --recipe melanoma-markers --max-cells 400 --seed 21
.venv-software/bin/vtsc run-tissue outputs/melanoma10x_input.zip outputs/melanoma_tissue --scenario declared-recognition --hours 24 --seed 21
.venv-software/bin/vtsc run-tissue outputs/melanoma10x_input.zip outputs/melanoma_recognition_off --scenario recognition-off --hours 24 --seed 21
```

Esses comandos antigos **não satisfazem** o novo objetivo anti-PD-1/estados.
Preserve-os para regressão ou migre com versionamento e mensagens claras.

### 5.5 Outras fontes já inspecionadas

- Sade-Feldman / GSE120575:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575 .
  `DATA/SadeFeldman2018/raw/Cell2018_TableS1.xlsx`, origem
  https://ars.els-cdn.com/content/image/1-s2.0-S0092867418313941-mmc1.xlsx ;
  `DATA/SadeFeldman2018/raw/GSE120575_patient_ID_single_cells.txt.gz`, origem
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz .
  Só esses arquivos foram adquiridos; não há matriz de expressão completa
  desse estudo pronta neste projeto. É Smart-seq2, CD45+, tempos posteriores.
  Resumos publicados de estados foram expostos em uma inspeção anterior:
  não chamar essa coorte de teste final intocado.
- `DATA/Liu2026/raw/biorxiv-2026.09.12.751186v1.pdf`, origem
  https://www.biorxiv.org/content/10.64898/2026.09.12.751186v1.full.pdf .
  Artigo baixado, não matriz celular; disponibilidade de dados não estabelecida.
- GSE272993 e GSE164237: sangue/PBMC/CD8 periférico; não substituem tumor.
  GSE123813: BCC/SCC, não melanoma. Há snapshots de registros, não seus datasets
  completos, em `reports/evidence/cohort-expansion/`.
- Fontes/snapshots, checksums e status de arquivos parciais estão em
  `data/COHORT_EXPANSION_MANIFEST.json`. Não trate `.part` ou resposta HTML como dados.

## 6. O que a auditoria da coorte já demonstrou

Não repita descobertas como se fossem novas e não confunda limite clínico com
pares disponíveis no objeto:

- SPECIAL: 3 candidatos em monoterapia com ≥30 CD8 em ambas as biópsias;
  6 com ≥15; 3 com ≥50. São limites por contagem, não elegibilidade definitiva.
- SPECIAL + pacientes novos de Lin: no limiar principal ≥30, **4 Mono e 8 Combo**.
  Total 12, abaixo do piso anterior de 10 desenvolvimento +4 teste; não é
  autorização para agrupar tratamentos. Os tetos clínicos antigos 9/23 foram
  superados pela auditoria do objeto.
- Lin tem 25 amostras internas/17 pacientes, apenas 8 com Pre e On depositados;
  a tabela clínica tinha mais pares. P2 tem duas datas On possíveis sem vínculo
  resolvido com a biblioteca. Há reutilização de SPECIAL no atlas.
- SPECIAL: 26/46 biópsias eram de linfonodo; local é covariável importante.
  Mesma lesão não foi verificada. Química mista tem divergência entre fontes
  EGA e S1. Não fixar uma interpretação conveniente.
- As duas S1 **não forneceram os intervalos individuais esperados**. Descrições
  gerais de “precoce”/ciclo 2 não autorizam inventar Δt por paciente.
- Lockbox não foi atribuído, não houve ajuste de modelo nem análise nova de
  mudanças individuais de estado SPECIAL. Quatro pacientes não permitem
  rejeição Wilcoxon exata unilateral a 0,05 (mínimo 0,0625).
- TME completo pode revelar células adicionais, mas isso precisa ser demonstrado
  por regra de linhagem/QC fixada antes da contagem, sem otimização para passar G1.

Na entrega exploratória, use o 10x independente para desenvolver comportamento
simulado. Não consuma a futura avaliação SPECIAL para escolher regras e efeitos.

## 7. Ambientes, reprodução e referência Xenium

No checkout local:

- `.venv-software/bin/python`: ambiente do produto/notebooks; Python 3.12.13.
  `pyproject.toml` aceita 3.12/3.13; dependências principais pinadas.
- `.venv/bin/python`: ambiente científico anterior, `requirements.lock.txt`,
  `environment.yml`. Interface CmdStanPy existe; backend de inferência ainda
  não foi decidido e compilador/modelo Stan não foi instalado para este estudo.
- R 4.6.0; `.Rlib`; `renv.lock`, `sessionInfo.txt`. Leitura/exportação usou
  Seurat 5.6.0, SeuratObject 5.4.0 e Matrix 1.7.6.
- Extensões Lin: `.Rlib-expansion`, `.pydeps-expansion`,
  `renv-expansion.lock`, `requirements-expansion.lock.txt`.
- Restauração: `r/restore_environment.R`, `r/restore_expansion_environment.R`;
  `r/setup_environment.R` documenta preparação anterior.

Primeira verificação, sem reprocessar pacientes:

```bash
cd /Users/heldernakaya/Desktop/2026/Papers/VirtualTissue/scRNAseq/VirtualTissue-SC
.venv-software/bin/vtsc verify-freeze
MPLCONFIGDIR=/tmp/vtsc-matplotlib .venv-software/bin/python -m pytest -q
```

Em uma instalação limpa, use Python compatível e `pip install -e '.[test,notebooks]'`.
Para auditoria clínica, o extra `[research]` acrescenta leitores PDF/XLSX, mas
não substitui o ambiente científico completo. Registre mudanças de dependências
sem sobrescrever os locks históricos. Não presumir que testes macOS provam
funcionamento de Windows/Linux/Colab. Kernels podem precisar de permissões para
portas locais; use a infraestrutura disponível, sem diagnosticar isso como falha biológica.

Arquivos reutilizáveis do produto: `src/vtsc/{inputs,model,tissue_inputs,tissue,
tissue_results,results,provenance,cli}.py`, `src/vtsc/web/tissue.html`, `tests/`,
`notebooks/Part1_Prepare_scRNAseq.ipynb`, `Part2-3_Run_Simulation.ipynb`,
`Part4_Explore_Results.ipynb`, `scripts/validate_notebooks.py`.

Referência XENIUM, somente leitura:

```text
XENIUM/src/virtualtissue_xenium/engine/xenium_runtime/v2/
  mechanics.py, menus.py, local_policy.py, autonomous.py, fields.py,
  scheduler.py, provider_jev.py, geometry.py, runtime.py
XENIUM/src/virtualtissue_xenium/engine/web/player_template.html
XENIUM/src/virtualtissue_xenium/engine/player.py
```

Reaproveitar conceitos/trechos exige manter autoria/licença quando aplicável.
Não transportar como medidas as geometrias Xenium para dados dissociados.

## 8. Contrato mínimo do modelo que você deve construir

Antes de implementar, escreva uma especificação curta com equações, unidades,
proveniência, parâmetros e limites. Consulte fontes primárias para sustentar
mecanismos; cite o que realmente suporta cada hipótese, sem inventar taxas.

1. **Identidade e estado separados.** Cada agente tem `cell_id`, `patient_id`/
   `sample_id` quando conhecidos, `lineage`, `cd8_state`, `time_in_state`,
   `parent_id`, `clone_id` se conhecido/simulado, `avatar/reference_id` e incerteza.
   Clone TCR não implica especificidade antigênica ou árvore mãe-filha observada.
2. **Estados e eixos separados.** Comece pelas quatro categorias publicadas
   Memory, Activated, Cytotoxic, Dysfunctional, auditando cobertura e separabilidade.
   Represente disfunção/exaustão, proliferação, citotoxicidade e persistência como
   eixos distintos. Uma célula disfuncional pode proliferar. Não imponha uma cadeia
   única irreversível nem trate TCF7 isolado como diagnóstico de Tpex.
3. **Grafo explícito.** Defina transições permitidas e proveniência de cada aresta.
   No estudo, usar robustez/treinamento e dados pareados conforme plano. No demo,
   grafo assumido deve estar rotulado; não alegar que foi inferido da única amostra.
4. **Molecular observado versus simulado.** Preserve RNA inicial imutável.
   Mudanças de estado podem selecionar um perfil de referência/alterar capacidades
   simuladas, mas não criar “expressão medida futura”. Mostre as duas camadas.
5. **Núcleo de hazards comum.** O mesmo modelo de estados deve alimentar o solver
   não espacial e os agentes espaciais. Documente probabilidades válidas por passo,
   transições concorrentes, crescimento, morte e influxo. Evite dois motores
   independentes com comportamento ajustado separadamente.
6. **Modelo não espacial.** Implemente M0 persistência, M1/M1b mudança média,
   M2 compartimentos com transição/crescimento líquido/influxo e M3 condicionado
   a ambiente. Implemente verossimilhança/estimação e recuperação sintética,
   mesmo quando o ajuste real permaneça indisponível. Parâmetros assumidos e
   estimados precisam ter tipos/proveniência e saídas diferentes.
7. **Anti-PD-1 explícito.** Defina uma exposição temporal/schedule configurável
   e um termo de checkpoint/reconhecimento que modifique mecanismos/hazards
   especificados. Separe presença de alvo, reatividade, contato, sinalização e
   capacidade efetora. Expressão de PDCD1/CD274 não é medida automática de
   proteína, ocupação ou resposta. Não impor que anti-PD-1 reverte toda exaustão,
   transforma todas as células em efetoras ou sempre mata mais tumor.
8. **Unidades e incerteza.** Se não houver farmacocinética validada, use exposição
   normalizada claramente assumida; não apresentar dose clínica calibrada.
   Cada taxa/constante/fração de reatividade deve dizer se é estimada, inferida,
   assumida ou medida. Faça análise de sensibilidade dos parâmetros dominantes.
9. **Controles adequados.** Usar mesmos inputs/posições e seeds emparelhados para
   anti-PD-1 versus exposição zero, reconhecimento off, checkpoint desacoplado
   quando apropriado, transições off/frozen e contexto espacial embaralhado.
   Exposição zero é um **contrafactual simulado**, não uma biópsia controle observada.
10. **Espaço assumido.** Layout, densidade e escala devem ser escolhas registradas.
    Preserve composição; se houver enriquecimento artificial, declare-o. Implemente
    contato, movimento/exclusão, divisão, morte/clearance e exposições locais
    suficientes para a hipótese; não adicionar dezenas de campos sem evidência.
11. **Decisão rastreável.** Registre estado anterior/novo, opções permitidas,
    hazards/probabilidades, draw quando estocástico, tempo-no-estado, exposições,
    anti-PD-1, alvo/contatos, regra/provedor, resultado e seed. Uma explicação
    textual posterior não substitui o registro causal do motor.
12. **Jev preparado, não fingido.** Interface opcional para escolher transições
    permitidas, com mock/local policy e journal verificável. Não chamar a regra
    local de IA/Jev, não usar confiança do LLM como taxa física. Comparação
    científica e chamadas pagas continuam condicionadas aos critérios do plano.

## 9. Plano operacional: tarefas, dependências e aceitação

Crie um quadro versionado para as tarefas P abaixo. Preserve o quadro T/G do
estudo original. Estes são entregáveis a implementar, não coisas já feitas.
Cada tarefa precisa de arquivo de evidência, teste/inspeção pertinente, status
e critério de aceite. Não marcar fases inteiras como completas por uma captura de tela.

### Fase P0 — Retomada, contrato e segurança do histórico

- **P0.1 — Inventário e baseline.** Ler arquivos da seção 3, confirmar caminhos,
  conferir freeze e testes. Saída: relatório de baseline, comandos e inventário
  atualizado. Aceite: distinguir cada arquivo ausente, parcial, bruto, derivado
  ou quarentenado; preservar alterações preexistentes.
- **P0.2 — Emenda de escopo.** Depende P0.1. Registrar produto exploratório
  autorizado, limites do estudo, dados/saídas permitidos e proibição de contato
  externo. Aceite: não alterar freeze; sem alegar que gates científicos passaram.
- **P0.3 — Especificação e critérios.** Depende P0.2. Escrever contrato da seção 8,
  schema de input/run, backend, estado/eixos, regras de tratamento, controles,
  orçamento de memória/runtime e mapa de rastreabilidade. Aceite: cada efeito
  implementável e testável; cada parâmetro com unidade/proveniência.

### Fase P1 — Dados e entrada do usuário

- **P1.1 — Exemplo melanoma e anotação.** Depende P0.3. Reusar matriz 10x verificada,
  auditar QC/cobertura e classificador de linhagem/estados. Aceite: relatório com
  incerteza, Unknown e composição; sem adição silenciosa de células.
- **P1.2 — Adapter SPECIAL/TME.** Depende P0.2. Fixar auditoria de linhagem/QC
  sem desfecho; adquirir os dois RDS públicos faltantes, verificar, inspecionar,
  exportar e validar round-trip em diretórios novos. Aceite: matrizes alinhadas
  por barcode, manifest e restrições de acesso. Não ajustar regras para passar G1.
- **P1.3 — Lin e fontes complementares.** Depende P0.2. Manter status metadata-only
  até validar expressão/alinhamento; resolver o que for possível publicamente,
  excluir reutilização de SPECIAL dos novos pacientes. Aceite: nenhuma matrícula
  clínica vira par celular por suposição. Se a memória impedir a leitura, registrar
  a pendência e continuar produto/10x; não repetir conversão de 29 GB inutilmente.
- **P1.4 — Own-data completo.** Depende P1.1 e schema P0.3. H5AD/raw layers,
  Cell Ranger H5 e MTX, metadados por barcode, uma biópsia por execução,
  exportação Seurat documentada, annotations revisadas ou receita explícita.
  Aceite: contagens/orientação corretas, IDs únicos, erros claros para camada
  normalizada, amostras misturadas, gene IDs sem mapping e linhagens/estados inválidos.

### Fase P2 — Referência CD8, programas e grafos

- **P2.1 — Biblioteca de programas.** Depende P1.1. Fontes/versões/cobertura,
  eixos independentes, sobreposição e controles de dissociação/profundidade.
  Aceite: biblioteca congelada; AUCell/UCell/alternativa justificada; não
  reutilizar os três escores ilustrativos como substituto sem avaliação.
- **P2.2 — Estados e incerteza.** Depende P2.1. Mapear quatro estados CD8 com
  confiança/Unknown; separar cycling e contexto. Aceite: validação com marcadores
  reservados e procedência do reference; nada ajustado nos pacientes de teste.
- **P2.3 — Grafo e perfis de estado.** Depende P2.2. Arestas, direção quando
  defensável, tabela de capacidades e perfis simulados. Aceite: grafo assumido
  e grafo inferido identificados; nenhum pseudotempo apresentado como tempo real.
- **P2.4 — Ambiente e covariáveis.** Depende P1.4/P2.1. Features tumorais/mieloides,
  local, química e denominadores corretos. Aceite: dados imunes isolados não
  produzem “fração de tumor no TME”; features ausentes não são zeros inventados.

### Fase P3 — Dinâmica e anti-PD-1

- **P3.1 — Solver comum de estados.** Depende P2.3/P2.4. Generator/hazards,
  crescimento/influxo e observação de composição. Aceite: conservação quando
  aplicável, probabilidades válidas, casos analíticos e testes de extremos.
- **P3.2 — Ladder e estimador.** Depende P3.1. M0/M1/M1b/M2/M3, modelo de
  observação Dirichlet-multinomial, incerteza de baseline, regularização/priors.
  Aceite: recuperação em verdade sintética e relatórios de identificabilidade;
  backend pinado. Ajuste real depende de elegibilidade, não de vontade de mostrar ganho.
- **P3.3 — Mecanismo anti-PD-1.** Depende P3.1 e contrato P0.3. Implementar
  schedule/exposição/checkpoint e como afetam cada hazard/ação. Aceite: exposição
  zero reproduz o baseline correspondente; desativar mecanismo desacopla a droga;
  ausência de alvo/reatividade obedece às regras declaradas; sem benefício garantido.
- **P3.4 — Estados funcionais e histórias.** Depende P3.3. Tempo-no-estado,
  diferenciação, disfunção/exaustão e efeitos nas capacidades/ações. Aceite:
  transições alteram o motor, não só cor/legenda; fadiga, disfunção e proliferação
  continuam distintas; morte irreversível e herança na divisão são especificadas.
- **P3.5 — Cenários e ensembles.** Depende P3.2–P3.4. Gerar arms emparelhados,
  controles e sensibilidade com seeds. Aceite: incerteza/variabilidade visível,
  não escolher o seed com efeito mais bonito; contrafactual rotulado.

### Fase P4 — Tecido espacial derivado do mesmo modelo

- **P4.1 — Runtime de estado dinâmico versionado.** Depende P3.4. Evoluir schema
  e agentes sem quebrar silenciosamente bundles antigos. Aceite: linhagem,
  estado, clone, parent, avatar e fonte molecular separados; referências preservadas.
- **P4.2 — Mecânica e exposição local.** Depende P4.1. Reusar/refinar movimento,
  contato, crescimento/divisão e campo mínimo. Aceite: limites espaciais/velocidade,
  resolução temporal, localização do tratamento/campo e contabilidade corretos.
- **P4.3 — Consistência com não espacial.** Depende P4.2. Fixar inputs locais
  nas médias de amostra e comparar ensembles ao solver comum. Aceite: tolerância
  justificada, nenhuma equivalência baseada numa única seed, diferenças espaciais
  quantificadas; incluir mecanismos extras ou desativá-los coerentemente na comparação.
- **P4.4 — Geometrias e controles.** Depende P4.3. Layouts relevantes/assumidos,
  inicialização proporcional estratificada explícita e embaralhamento. Aceite:
  distinguir número de células, composição, geometria e seed; não chamar posição
  assumida de reconstrução histológica.
- **P4.5 — Interface de política.** Depende P4.1. Estocástica/hazard, política
  local e mock de provedor; menus legais e logs. Aceite: nenhum provedor executa
  ações fora do contrato; offline é completo; chamadas Jev reais não são requisito
  para terminar a entrega atual.

### Fase P5 — Experiência de uso e inspeção

- **P5.1 — Viewer de células/estados.** Depende P4.2. Células distinguíveis,
  movimentos/contactos verdadeiros do registro, fronteiras/núcleos esquemáticos,
  cores por linhagem OU estado, histórico de divisão/morte/transição. Aceite:
  renderer não inventa eventos e não modifica o motor; geometrias assumidas visíveis.
- **P5.2 — Inspetor de decisões.** Depende P5.1/P3.4. Clique/seleção acessível,
  estado anterior/novo, tempo-no-estado, quatro eixos, opções, probabilidades,
  sorteio, contexto, drug exposure, checkpoints e provenance. Aceite: cada item
  corresponde ao journal; distinguir RNA t0 observado e perfil de estado simulado.
- **P5.3 — Comparação anti-PD-1/controle.** Depende P3.5/P5.2. Reprodução/scrub
  sincronizados, composição CD8, transições, proliferação, morte, tumor e exposição.
  Aceite: braços mesmo input/contrato/condições iniciais; diferenças de config
  explícitas; curvas e contagens consistentes com os arquivos.
- **P5.4 — Entrada → execução → resultados.** Depende P1.4/P5.3. CLI e três
  notebooks, config de cenário, retomada/export ZIP e viewer offline. Aceite:
  usuário executa dados próprios sem editar código do motor; exemplos pequenos
  rápidos; reruns não sobrescrevem silenciosamente resultados aceitos.

### Fase P6 — Verificação, documentação e pacote GitHub

- **P6.1 — Testes científicos de implementação.** Depende P3/P4. Testar limites,
  transições/hazards, tratamento, controles negativos, recuperação/identificabilidade,
  equivalência média, morte/divisão/influxo, leakage guards, integridade e sementes.
  Aceite: testes verificam consequências, não só repetem fórmulas de implementação.
- **P6.2 — Integração e inspeção visual.** Depende P5/P6.1. Rodar do zero o
  fixture sintético, exemplo real melanoma e fixture own-data; executar notebooks,
  reabrir ZIP e testar UI no navegador. Aceite: evidências/screenshots/logs,
  erro útil para entrada inválida e teste real de wheel com HTML incluído.
- **P6.3 — Documentação final.** Depende P6.2. README com instalação/quickstart,
  own-data, mecanismo anti-PD-1, estados, limites, troubleshooting, outputs,
  licenças/attribution/citação; ODD; mapas tarefa→evidência e afirmação→arquivo.
  Aceite: nada promete clínica validada ou links de execução ainda não testados.
- **P6.4 — Release local testável.** Depende P6.3. Versão nova coerente, locks,
  changelog, CI útil, pacote pequeno, exemplos/downloader, ignores de dados e
  credenciais, instruções de clone/push. Aceite: checkout limpo instalável e
  auditado; diff revisado; commits locais; sem push ou publicação automática.

## 10. O plano científico original: 46 entradas preservadas

A lista abaixo é o mapa de continuidade, não autorização para declarar tarefas
pendentes como prontas após construir o produto. O plano original contém
critérios detalhados, dependências e suplemento associado; leia a versão
congelada para executá-los. Contagem do quadro em 28/10: 6 concluídas,
2 parciais, 38 pendentes. “Concluída” é o status registrado da auditoria anterior.

| ID | Tarefa/gate do plano original | Status registrado |
|---|---|---|
| T0.1 | Repository and environment | Concluída |
| T0.2 | Data acquisition and integrity | Concluída |
| T0.3 | Licence and data-use note | Concluída |
| T0.4 | Pre-registration freeze | Concluída |
| T1.1 | Export Seurat objects to a neutral format | Concluída |
| T1.2 | Patient and sample table | Concluída |
| T1.3 | Eligibility and flow diagram | Parcial |
| T1.4 | Lockbox assignment | Pendente |
| T1.5 | CD8 state availability (development patients only) | Pendente |
| T1.6 | Technical confounders | Parcial |
| T1.7 | Feasibility report · Gate G1 | Pendente |
| T2.1 | CD8 subset QC and annotation audit | Pendente |
| T2.2 | Program library and coverage | Pendente |
| T2.3 | Scoring and fold-aware calibration | Pendente |
| T2.4 | State graph, pseudotime and robustness | Pendente |
| T2.5 | Functional axes per state | Pendente |
| T2.6 | Environment features | Pendente |
| G2 | State reference frozen | Pendente |
| T3.1 | Implement the model ladder | Pendente |
| T3.2 | Simulation-based parameter recovery | Pendente |
| T3.3 | Identifiability on real data | Pendente |
| T3.4 | Development LOPO evaluation | Pendente |
| T3.5 | Sensitivity analyses | Pendente |
| T3.6 | Freeze the model · Gate G3 | Pendente |
| T4.1 | Lockbox evaluation | Pendente |
| T4.2 | External direction check (GSE120575) | Pendente |
| T4.3 | Complementary temporal and clonal data (optional) | Pendente |
| G4 | Go/no-go for spatial work | Pendente |
| T5.1 | SPECIAL Xenium audit | Pendente |
| T5.2 | Runtime design document (ODD) and margins | Pendente |
| T5.3 | Implement the dynamic-state runtime | Pendente |
| T5.4 | Mean-field consistency test | Pendente |
| T5.5 | Geometry experiments | Pendente |
| T5.6 | Spatial evaluation · Gate G5 | Pendente |
| T6.1 | Pilot (≤ 100 calls, explicit approval) | Pendente |
| T6.2 | Policy comparison · Gate G6 | Pendente |
| T7.1 | Combination-therapy pairs (separate model or test of transfer) | Pendente |
| T7.2 | Site-specific templates (lymph node vs skin) beyond the covariate | Pendente |
| T7.3 | Healthy-skin context (needs a tissue profile, new data, new menus) | Pendente |
| T7.4 | Explicit untreated counterfactual arm (labelled hypothetical) | Pendente |
| T7.5 | Coupling to VirtualTissue-LN (priming compartment → influx) | Pendente |
| T8.1 | Compile supplementary notes | Pendente |
| T8.2 | Figures | Pendente |
| T8.3 | Data and code availability | Pendente |
| T8.4 | Reproducibility run | Pendente |
| T8.5 | Claims audit | Pendente |

Restrições científicas que permanecem:

- G1: ≥10 pacientes de desenvolvimento elegíveis e disponibilidade dos estados;
  manter split por paciente. Não usar células como replicações independentes.
- T2: referência/grafo/calibração apenas no treino de cada fold; não usar
  embeddings integrados dos autores como representação ajustada independente.
- T3: recuperar parâmetros/avaliar identificabilidade antes de interpretá-los;
  comparar modelos simples, usar métricas composicionais e incerteza.
- G3: congelar modelo antes do teste; G4 controla inferência científica sobre
  valor do espaço; G5 exige desempenho e padrão espacial predefinido.
- G6: Jev exige comparação, rótulos cegos/controles e orçamento aprovado.
- Monoterapia inclui nivolumabe e pembrolizumabe. Combinação tem estrato próprio;
  não assumir efeitos iguais nem contraste causal entre tratamentos não randomizados.
- Não prever resposta clínica; desfechos futuros não entram no baseline.
- Resultados exploratórios/piloto e resultados confirmatórios não compartilham
  rótulos de evidência. Um modelo estimado com n pequeno não vira validado por rodar.

Suplemento a manter alinhado: **S1** fontes/licenças; **S2** inventário;
**S3** coorte/confundidores; **S4** referência/escores/eixos; **S5** grafo;
**S6** ambiente/modelos; **S7** splits/inferência; **S8** recuperação/identificabilidade;
**S9** sensibilidades; **S10** ODD/consistência espacial; **S11** validação externa;
**S12** software/Jev/custo; **S13** limites das afirmações. Preserve
`tables/claims_ledger.tsv` e gere figuras somente de arquivos de resultado versionados.

## 11. O que entregar ao encerrar esta implementação

Só considere a entrega de engenharia completa quando houver:

1. Exemplo real de melanoma 10x inicializando tumor e microambiente; fixture
   sintética pequena e entrada own-data passando ponta a ponta.
2. CD8 mudando entre estados permitidos com registro de probabilidades e
   capacidades atualizadas, inclusive disfunção/exaustão distinta de fadiga.
3. Anti-PD-1 como mecanismo/schedule explícito, controles coerentes e comparação
   entre braços; exposição zero não confundida com reconhecimento off.
4. Motor não espacial e espacial ligados pelo mesmo núcleo, com consistência
   verificada em condições comparáveis e limitações documentadas.
5. Viewer com células interagindo e transições reais, timeline, inspeção individual,
   campo/tratamento, comparação de arms e separação observado/inferido/assumido.
6. Logs reprodutíveis, resultados portáveis, verificações relevantes e screenshots.
7. README, três notebooks, ODD, licenças, schemas, exemplos, CI e release local
   prontos para o grupo instalar e usar; sem dados grandes/privados no Git.
8. Tabela honesta do que continua faltando para a validação SPECIAL original,
   sem impedir que eu teste o software exploratório entregue.

Comece agora lendo os arquivos e verificando o baseline. Em seguida registre o
contrato/progresso e implemente P0–P6. Não repita o desvio PBMC, não pare depois
de melhorar apenas o visual e não me devolva outro plano em lugar do software.
