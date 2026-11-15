
# Burrows–Wheeler Transform (BWT) & FM-Index:

If you’ve ever wondered how bioinformatics tools like **Bowtie2** or **BWA** can search through the entire human genome (3.2 billion DNA letters!) on a regular laptop without it catching fire, the answer is the **Burrows–Wheeler Transform (BWT)**.

Normally, if you want to find a short DNA read inside a massive genome, scanning from start to finish takes $O(N)$ time—which is way too slow when you have millions of reads. BWT solves this by scrambling a sequence in a very specific, **100% reversible** way so that:
1. **Identical letters group together into "runs"** (like `AAAA` or `TTTT`), making the genome super easy to compress.
2. **Repeated patterns line up alphabetically**, allowing us to search for any pattern of length $M$ in just **$O(M)$ steps**—no matter how huge the genome is!

---

## How the Code Works (Step-by-Step)

### 1. Input Validation & Auto-Detection (`validate_raw_sequence` & `validate_bwt_input`)
Before doing any math, the code acts like a bouncer at the door to make sure the input won't break the algorithm:

* **Auto-Detecting the Sequence Type:**
  * If the string only uses `A, T, G, C, N`, it’s flagged as **DNA**.
  * If it has `U` instead of `T` (`A, U, G, C, N`), it’s **RNA**.
  * If it uses standard amino acid letters (like `M, K, W, V, F`), it’s **PROTEIN**.
  * Anything else (like `"BANANA"`—which is legally required in every BWT assignment ever written) is treated as **TEXT**.
* **Protecting the Sentinel (`$`):**
  * BWT needs a special "end-of-string" marker called the **sentinel (`$`)**. We must make sure the user didn't sneak a `$` into their raw sequence.
* **Fixing the ASCII Betrayal (`_sentinel_sort_key`):**
  * By rule, `$` **must** be alphabetically smaller than every other character. In Python’s standard ASCII table, a space (`' '`, ASCII 32) is actually smaller than `'$'` (ASCII 36). Our helper function `_sentinel_sort_key` forces `$` to have priority `0` so it *always* sorts to the very top row, no matter what characters are in your string.

---

### 2. Forward BWT — Scrambling the Sequence (`bwt_encode`)
How do we encode a sequence like `GATTACA` or `BANANA`? Let's trace `BANANA`:

1. **Slap the Sentinel (`$`) on the end:**  
   `BANANA` $\rightarrow$ `BANANA$` (Total length $N = 7$).
2. **Generate all $N$ Cyclic Rotations:**  
   We pluck letters off the front and move them to the back, one by one, until we’ve gone full circle. Notice how in every single row, **the Last letter (`L`) is actually the letter that came right BEFORE the First letter (`F`) in the original word!**
3. **Sort the Rotations Alphabetically (with `$` first):**  
   Sorting all the rotations groups identical suffixes ("right-contexts") together:

   | Row | Shift (`SA[i]`) | First (`F`) | Sorted Rotation Matrix | Last (`L` / BWT) | Note |
   | :-: | :-: | :-: | :--- | :-: | :--- |
   | **0** | 6 | `$` | `$BANANA` | **`A`** | |
   | **1** | 5 | `A` | `A$BANAN` | **`N`** | |
   | **2** | 3 | `A` | `ANA$BAN` | **`N`** | |
   | **3** | 1 | `A` | `ANANA$B` | **`B`** | |
   | **4** | 0 | `B` | `BANANA$` | **`$`** | **$\leftarrow$ Primary Index (Row 4)** |
   | **5** | 4 | `N` | `NA$BANA` | **`A`** | |
   | **6** | 2 | `N` | `NANA$BA` | **`A`** | |

4. **Extract the Outputs:**
   * **BWT String (`L` column):** Reading the last column top-to-bottom gives **`ANNB$AA`**. Look at how the `N`s and `A`s clumped together!
   * **Primary Index (`4`):** The row where our original, unshifted string (`BANANA$`, shift 0) ended up. Because the original string ends with `$`, `L[primary_index]` is *always* where `$` sits!
   * **Suffix Array (`SA = [6, 5, 3, 1, 0, 4, 2]`):** Just a list of how many characters each sorted row was shifted from the original string. We save this for genome coordinate lookups later.

---

### 3. Inverse BWT (`bwt_decode` & `build_lf_mapping_structures`)
Here is the wildest part of BWT: if someone deletes the entire rotation matrix and *only* hands you the BWT string **`ANNB$AA`** and the Primary Index **`4`**, how do you get `BANANA` back without rebuilding the $N \times N$ grid?

We use three clever pieces of information:

1. **Recovering the First Column (`F`):**  
   Since `F` is just all the characters of the sequence sorted alphabetically, we get `F` simply by sorting `L` (`ANNB$AA` sorted $\rightarrow$ `$AAABNN`).
2. **The C-Table (`C[char]`):**  
   `C[char]` counts how many total characters in the string are strictly smaller than `char`.  
   *Why do we care?* Because in a 0-based table, that count is the **exact row number where `char` starts in column `F`!**
   * For `BANANA$`: `$` starts at row `0` (`C['$'] = 0`), `A` starts at row `1` (`C['A'] = 1`), `B` starts at row `4` (`C['B'] = 4`), and `N` starts at row `5` (`C['N'] = 5`).
3. **The LF-Mapping Property (Last-to-First):**  
   In both `F` and `L`, multiple `'A'`s appear. How do we know *which* `'A'` in `L` matches *which* `'A'` in `F`?  
   * **Rule:** The **1st `'A'`** in `L` is the exact same text character as the **1st `'A'`** in `F`. The **2nd `'A'`** in `L` is the **2nd `'A'`** in `F`, and so on!  
   * *Why?* Because all rows starting with `'A'` in `F` are sorted by the letters that come *after* `'A'`. And all rows ending with `'A'` in `L` are *also* sorted by the exact same letters coming *after* `'A'`! Since both groups use the exact same tie-breaker, their relative order (`rank`) stays identical.

#### The Reconstruction Formula:
$$\text{LF}[i] = C[\,L[i]\,] + \text{rank}_0(L[i])$$
*(Translation: To find where the letter at `L[i]` lives in column `F`, take the starting row of that letter in `F` (`C[L[i]]`) and add how many times you've already seen that letter in `L` above row `i`.)*

* **Walking Backwards in $O(N)$ Steps:**
  1. We start at Row `0` (where `F[0] = '$'`, which is the end of our word).
  2. The letter in `L[0]` is `'A'` ( specifically the 1st `'A'` in `L`, so `rank = 0`). That means `'A'` comes right before `'$'`!
  3. We jump to `LF[0] = C['A'] + 0 = 1 + 0 = Row 1`.
  4. At Row `1`, `L[1]` is `'N'` (1st `'N'` in `L`, `rank = 0`), so `'N'` comes before `'A'`. We jump to `LF[1] = C['N'] + 0 = 5`.
  5. Repeat $N$ times until we hit `'$'` again! That spells `ANANAB` backwards, which we flip to get **`BANANA`**.

---

### 4. Biological Analysis & Repeat Clustering (`analyze_biological_sequence`)
Why do biologists care so much about BWT? This function bridges the computer science with real-world genomics:

* **GC / AT / AU Content:**
  * In DNA and RNA, **Guanine (`G`)** and **Cytosine (`C`)** form **3 hydrogen bonds** with each other, while **Adenine (`A`)** and **Thymine (`T`) / Uracil (`U`)** only form **2**.
  * A higher **GC%**\ means the DNA is much more stable and requires a higher temperature ($T_m$) to melt apart during lab techniques like PCR.
* **Homopolymers & Repeated $k$-mers (`find_repetitive_patterns`):**
  * Genomes are full of repeats—both **homopolymers** (like `AAAAA`) and **repeated $k$-mers** (short words of length 2, 3, or 4 like `ATG` that appear multiple times).
* **Run-Length Encoding (`count_runs`):**
  * When a motif like `ANA` repeats in a sequence, every rotation starting with `NA...` gets sorted into adjacent rows in the BWT matrix. That forces the letters right before `NA` (the `'A'`s) to sit next to each other in the BWT output!
  * `count_runs` compresses consecutive identical letters into counts (e.g., `BANANA` has **6 runs**: `B1 A1 N1 A1 N1 A1`, whereas its BWT `ANNB$AA` has fewer runs because the `N`s and `A`s grouped up: `A1 N2 B1 $1 A2`). On a real genome, this shrinks file sizes dramatically.

---

### 5. FM-Index Backward Search (`fm_index_backward_search`)
Instead of decoding the whole genome just to find where a short DNA read (pattern $P$, like `"ANA"`) is located, the **FM-Index** searches the compressed BWT directly by reading the pattern **backwards** (right-to-left):

1. **Keep an Occurrence Checkpoint Table (`Occ(char, row)`):**  
   This table simply records how many times `char` has appeared in the BWT column `L` from row `0` down to `row`.
2. **Narrow the Row Range `[start_row .. end_row]`:**  
   We start with all rows (`0` to $N-1$) as candidates, and read our pattern `"ANA"` from right to left (`'A'`, then `'N'`, then `'A'`):
   * At each letter `c`, we update our row range using the exact same LF-mapping logic:

$$
\mathrm{new\_start} = C[c] + \mathrm{Occ}(c, \mathrm{start\_row} - 1)
$$

$$
\mathrm{new\_end} = C[c] + \mathrm{Occ}(c, \mathrm{end\_row}) - 1
$$
   * **Step 1 (`'A'`):** Narrows `[0..6]` down to rows `[1..3]` (all rows starting with `'A'`).
   * **Step 2 (`'N'`):** Checks how many `'N'`s appear in `L[1..3]` (there are two!) and jumps to their rows in `F`: `[5..6]` (rows starting with `"NA"`).
   * **Step 3 (`'A'`):** Checks how many `'A'`s appear in `L[5..6]` (two!) and jumps to their rows in `F`: **`[2..3]`** (rows starting with `"ANA"`).
3. **Read the Genomic Coordinates from the Suffix Array (`SA`):**  
   Our final row range is `[2..3]`. Looking up `SA[2]` and `SA[3]` gives **`[3, 1]`**—meaning `"ANA"` appears at **index 1** and **index 3** in `"BANANA"`! And we found both matches in just 3 tiny math steps.

---              
## Implementation & AI Usage Note

### Implementation Note
* The core algorithmic workflows in this project—including cyclic rotation sorting, $C$-table construction, Last-to-First (LF) mapping inversion, and FM-Index backward search—were manually designed and implemented using the primary literature and educational coding resources listed in the references (specifically Burrows & Wheeler, 1994; Ferragina & Manzini, 2000; Compeau & Pevzner, 2018; the Rosalind `BA9I–BA9N` problem track; Prof. Ben Langmead’s JHU Python notebooks; and GeeksforGeeks).
* AI assistance was used strictly for post-implementation code tuning and debugging—such as fixing edge-case input validation, handling the custom ASCII sort priority for the `'$'` sentinel, and tidying up the CLI table formatting. All core algorithm design, biological sequence analyses, and test-case verifications were completed manually.

---

## References

### 1. Core Algorithms (BWT, LF-Mapping, Suffix Arrays, & FM-Index)
1. **Burrows, M., & Wheeler, D. J. (1994).** *A block-sorting lossless data compression algorithm.* Digital Equipment Corporation, SRC Research Report 124.  
   *(Original paper introducing the Burrows–Wheeler Transform, cyclic rotations, and the Last-to-First [LF] mapping property used in `bwt_encode` and `bwt_decode`.)*
2. **Ferragina, P., & Manzini, G. (2000).** *Opportunistic data structures with applications.* In Proceedings of the 41st Annual Symposium on Foundations of Computer Science (FOCS) (pp. 390–398). IEEE. https://doi.org/10.1109/SFCS.2000.892127  
   *(Introduced the **FM-Index**, combining the BWT with the $C$-table and Occurrence [$Occ$] checkpoints for $O(|P|)$ backward search used in `fm_index_backward_search`.)*
3. **Manber, U., & Myers, G. (1993).** *Suffix arrays: A new method for on-line string searches.* SIAM Journal on Computing, 22(5), 935–948. https://doi.org/10.1137/0222058  
   *(Introduced the Suffix Array [$SA$] data structure used to recover genomic coordinates from BWT row intervals.)*

### 2. Bioinformatics & Genome Alignment Applications
4. **Li, H., & Durbin, R. (2009).** *Fast and accurate short read alignment with Burrows–Wheeler transform.* Bioinformatics, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324  
   *(Paper behind **BWA**, demonstrating how BWT and FM-Index allow ultra-fast mapping of DNA sequencing reads to large reference genomes.)*
5. **Langmead, B., & Salzberg, S. L. (2012).** *Fast gapped-read alignment with Bowtie 2.* Nature Methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923  
   *(Describes **Bowtie 2**, combining FM-Index backward search with dynamic programming for gapped DNA read alignment.)*
6. **Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019).** *Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype.* Nature Biotechnology, 37(8), 907–915. https://doi.org/10.1038/s41587-019-0201-4  
   *(Describes **HISAT2**, which uses hierarchical BWT/FM-Index structures for fast RNA-Seq splice-junction alignment.)*

### 3. Textbooks & Documentation
7. **Compeau, P., & Pevzner, P. A. (2018).** *Bioinformatics Algorithms: An Active Learning Approach* (3rd ed., Chapter 9: "How Do We Locate Disease-Causing Mutations?"). Active Learning Publishers. ISBN: 978-0990374633.  
   *(Standard bioinformatics textbook covering BWT run-length clustering, the First-Last property, and pattern matching in genomes.)*

### 4. Web/code references
9. **Langmead, B.** *Computational Genomics Course Notebooks: Burrows–Wheeler Transform, Reversing BWT, and FM Index (Python).* Johns Hopkins University. Available at: https://github.com/BenLangmead/comp-genomics-class
10. **Rosalind Bioinformatics Platform.** *Bioinformatics Textbook Track (Problems BA9I–BA9N: BWT Construction, LF-Mapping, and BWMatching).* Available at: https://rosalind.info/problems/list-view/?location=bioinformatics-textbook-track
11. **GeeksforGeeks.** *Burrows–Wheeler Data Transform Algorithm & Inverting the BWT in Python.* Available at: https://www.geeksforgeeks.org/burrows-wheeler-data-transform-algorithm/
