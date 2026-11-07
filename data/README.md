# Data

Two CSV files, one row per provisioning visit, with nest-, group- and
watch-level values repeated on every row.

| File | Visits | Watches | Nests |
|---|---|---|---|
| `raw/mono_sheet_coordination.csv` | 1,561 | 113 | 63 socially monogamous |
| `raw/poli_sheet_coordination.csv` | 2,815 | 184 | 100 socially polyandrous |

Nests were filmed at the Dunedin Botanic Garden, New Zealand, when chicks were
8 and 9 days old. Visits during which a bird brooded, visits by unidentified
birds, and nests with unringed parents were excluded when the data were
compiled (see Methods of the manuscript).

## Columns used by the analysis

| Column | Meaning |
|---|---|
| `est.rep` | breeding season (1–7, 2009/10 to 2015/16) |
| `nest.id` | nest, within season |
| `sms` | social mating system (`mono` = pair, `polia` = one female and two males) |
| `female`, `male.a`, `male.b` | colour-ring combinations of the group's carers |
| `focal.id` | nest watch |
| `focal.date` | date of the watch |
| `focal.length.h` | effective watch duration, h (from the first visit to the end of the video) |
| `n.chicks` | number of chicks in the nest during the watch |
| `chick.age` | chick age, days |
| `visit.id` | visit number within the watch |
| `parental.status` | carer making the visit: `f` female, `a` alpha male, `b` beta male |
| `parental.sex` | `f` or `m` |
| `start.time`, `end.time` | arrival and departure, h from the start of the recording |

The remaining columns are per-watch summaries (visit counts, alternation counts
and rates, inter-visit intervals, synchrony) computed when the data were
compiled. The analysis does not use them; it recomputes every metric from the
visit sequences.

Three watches (`focal9`, `focal150`, `focal585`) have a visit out of
chronological order or a colour-ring combination that does not match the nest
roster; they are kept in the main analysis and dropped in a sensitivity
analysis.
