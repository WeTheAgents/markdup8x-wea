# Slack draft — markdup-wea parity handoff

> **Not sent yet.** Review, edit, paste into the channel when ready.
> Do not send until `docs/parity-evidence.md` and the README update
> have landed on `main` — the message links to them.

---

**Subject / thread:** markdup-wea — parity evidence ready, asking for one round of external testing

**Body:**

Команда, короткий апдейт по markdup-wea (Rust-реимплементация Picard MarkDuplicates, дроп-ин для `nf-core/rnaseq` шага `PICARD_MARKDUPLICATES`).

**Что сделано.**

1. **Non-UMI, default-off.** `samtools view | md5sum` байт-в-байт совпадает с Picard 3.4.0 на 8/8 ENCODE RNA-seq сэмплов: 1,549,871,416 записей, ~934M дубов, **0 flag-divergence**. 3.26× быстрее, 54× меньше peak RAM на K562_REP1 (325 MB vs 17.6 GB).
2. **UMI (`--barcode-tag RX`).** Валидирован на 3 независимых открытых датасетах, 2 организма, 3 геометрии:
   - GSE75823 (human UHRR, SCRB-seq 6nt UMI, 210M SE ридов)
   - PRJNA416930 (mouse testis, 8nt UMI-on-R2, 41.6M SE ридов)
   - GSE134031 (mouse microglia, QuantSeq inline 6nt UMI + TATA spacer, 3.78M SE ридов)
   Flag-set (сортированные `(qname, dup-bit)`): **0 lines diff на всех трёх**; metrics exact match. wea 2–2.5× быстрее Picard на двух замеренных.
3. **Regression gate.** После A.4-рефакторинга `SingleEndTracker` (окно-буфер для UMI-интерлива) перепрогнали 8-ENCODE — 0 divergence, 1.55B ридов. Любое изменение в `src/groups.rs` / `src/scan.rs` прогоняется через тот же gate.

**Чего нет в нашей доказательной базе.**

- Patterned flowcell (HiSeq X / NovaSeq) — единственный открытый корпус GTEx, dbGaP-gated (phs000424). Заявку **не подаём** — как индивидуалы в биоинформатике это социально дорогой путь.
- Single-cell (10x Chromium, Drop-seq) — out of scope, другая семантика дедупа.
- Multi-library BAMs — алгоритм поддерживает, парити-тестов нет.
- Long-read (PacBio / ONT) — не поддерживается; `FWD_WINDOW=1024` под Illumina, на более длинных ридах wea падает с понятным сообщением (fail-fast, не silent corruption).
- Optical duplicate detection — не реализовано (`deviations.md §1`).
- `DUPLEX_UMI` — не реализовано (`A.7`).

Полная доказательная база: [`docs/parity-evidence.md`](https://github.com/peachgabba22/markdup-wea/blob/main/docs/parity-evidence.md). Deviations: [`docs/deviations.md`](https://github.com/peachgabba22/markdup-wea/blob/main/docs/deviations.md).

**Вопросы / просьбы к вам.**

1. **Какие chemistries / protocols реально крутятся у вас / клиентов**, которых нет в базе выше? NovaSeq, multi-library, scRNA, targeted panel — что конкретно?
2. **Есть ли BAM (любой ваш, хоть и старый), на котором можно прогнать wea и прислать flag-diff?** Скрипт на 20 строк, на 50M ридов ~5 минут wall-clock. Готов приложить пошагово — нужно только согласие одного добровольца.
3. **Дальше куда двигаться?** Сейчас логические следующие треки — A.7 (DUPLEX_UMI + MC-tag parsing) и Phase 4 (multi-threaded Pass 1). Ок ли их пока отложить и переключиться на что-то более приоритетное для бизнеса до получения внешних подтверждений?

Спасибо за время. Любой ответ — даже «всё норм, продолжай» — уже полезен.

---

## Appendix — что именно попросить у добровольца

(Переводить в отдельное DM тому, кто согласится помочь с тестом.)

Нужен один coordinate-sorted BAM + соответствующий `*.markdup.bam` от
Picard 3.4.0 (либо мы прогоним Picard сами, если Picard-вывода нет).
Дальше всё автоматом:

```bash
# У нас на Hetzner:
WORK=/mnt/HC_Volume_105344878/tmp/external-<tag>
mkdir -p $WORK && cd $WORK
scp user@remote:/path/to/input.bam input.bam           # или rsync
scp user@remote:/path/to/picard.bam picard.bam         # если Picard у них уже прогнан

# Если Picard не прогнан, делаем сами — 4GB heap, ASSUME_SORTED=true:
java -Xmx4g -jar /opt/picard/picard.jar MarkDuplicates \
  I=input.bam O=picard.bam M=picard.metrics.txt \
  ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

# wea:
~/markdup-test/markdup-wea-a2/target/release/markdup-wea \
  --input input.bam --output wea.bam --metrics wea.metrics.txt
# добавить --barcode-tag RX если UMI

# diff:
samtools view picard.bam | awk '{print $1"\t"($2%2048>=1024)}' | sort > picard.flags
samtools view wea.bam    | awk '{print $1"\t"($2%2048>=1024)}' | sort > wea.flags
wc -l picard.flags wea.flags
diff picard.flags wea.flags | wc -l
diff picard.flags wea.flags | head -40     # если не ноль — это первые 40 расхождений
```

Ожидаемый результат: `diff ... | wc -l` == 0. Любой ненулевой остаток
→ мы берём разбор на себя, классифицируем, документируем.
