/// Output from the Hirschberg algorithm containing an optimal
/// global alignment and distance between two sequences.
#[derive(Debug, PartialEq, Eq)]
pub struct Output<'a, T, U> {
    alignment: Vec<(Option<&'a T>, Option<&'a U>)>,
    score: i32,
}

impl<'a, T, U> Output<'a, T, U> {
    /// Returns a slice of the optimal alignment. Each item is a tuple containing
    /// an element from sequence `a` or gap (indicated by `None`) at position 0, and
    /// an element from sequence `b` or gap at position 1. The ordering of the original
    /// sequences `a` and `b` is preserved.
    pub fn alignment(&self) -> &[(Option<&'a T>, Option<&'a U>)] {
        &self.alignment
    }

    /// Returns the score of the alignment.
    pub fn score(&self) -> i32 {
        self.score
    }

    fn combine(mut self, other: Self) -> Self {
        self.alignment.extend(other.alignment);
        Self {
            alignment: self.alignment,
            score: self.score + other.score,
        }
    }
}

/// Configuration for the Hirschberg algorithm
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_score: i32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_score: -1,
            gap_score: -1,
        }
    }
}

impl Config {
    /// Set the match score, which weights matching elements in the aligned sequences.
    pub fn match_score(mut self, score: i32) -> Self {
        self.match_score = score;
        self
    }

    /// Set the mismatch score, which weights differing elements in the aligned sequences.
    pub fn mismatch_score(mut self, score: i32) -> Self {
        self.mismatch_score = score;
        self
    }

    /// Set the gap score, which weights gaps in the alignment of sequences.
    pub fn gap_score(mut self, score: i32) -> Self {
        self.gap_score = score;
        self
    }

    pub fn compute<'a, T: PartialEq<U>, U: PartialEq<T>>(
        self,
        a: &'a [T],
        b: &'a [U],
    ) -> Output<'a, T, U> {
        if a.is_empty() {
            Output {
                alignment: std::iter::repeat(None)
                    .zip(b.iter().map(Option::Some))
                    .collect(),
                score: b.len() as i32 * self.gap_score,
            }
        } else if b.is_empty() {
            Output {
                alignment: a
                    .iter()
                    .map(Option::Some)
                    .zip(std::iter::repeat(None))
                    .collect(),
                score: a.len() as i32 * self.gap_score,
            }
        } else if a.len() == 1 || b.len() == 1 {
            self.needleman_wunsch(a, b)
        } else {
            let a_mid = a.len() / 2;
            let score_left = self.nw_score(&a[0..a_mid], b, false);

            let mut score_right = self.nw_score(&a[a_mid..], b, true);
            score_right.reverse();
            let b_mid = score_left
                .into_iter()
                .zip(score_right)
                .map(|(l, r)| l + r)
                .enumerate()
                .max_by_key(|(_, score)| *score)
                .unwrap()
                .0;

            let res = self
                .compute(&a[0..a_mid], &b[0..b_mid])
                .combine(self.compute(&a[a_mid..], &b[b_mid..]));
            res
        }
    }

    fn nw_score<T: PartialEq<U>, U: PartialEq<T>>(
        &self,
        a: &[T],
        b: &[U],
        reversed: bool,
    ) -> Vec<i32> {
        let mut row_1 = vec![0];
        let mut row_2 = vec![];

        for j in 0..b.len() {
            let next = row_1[j] + self.gap_score;
            row_1.push(next);
        }

        for i in 1..=a.len() {
            row_2.push(row_1[0] + self.gap_score);
            for j in 1..=b.len() {
                let sub = row_1[j - 1]
                    + self.replace(
                        &a[if reversed { a.len() - i } else { i - 1 }],
                        &b[if reversed { b.len() - j } else { j - 1 }],
                    );
                let del = row_1[j] + self.gap_score;
                let ins = row_2[j - 1] + self.gap_score;
                row_2.push(sub.max(del).max(ins));
            }
            row_1 = std::mem::take(&mut row_2);
        }
        row_1
    }

    fn replace<T: PartialEq<U>, U: PartialEq<T>>(&self, x: &T, y: &U) -> i32 {
        if x == y {
            self.match_score
        } else {
            self.mismatch_score
        }
    }

    fn needleman_wunsch<'a, T: PartialEq<U>, U: PartialEq<T>>(
        self,
        a: &'a [T],
        b: &'a [U],
    ) -> Output<'a, T, U> {
        let mut dp: Vec<Vec<Cell>> = (0..=b.len())
            .map(|row| {
                (0..=a.len())
                    .map(|col| Cell {
                        prev: if row == 0 && col == 0 {
                            None
                        } else {
                            Some((row.saturating_sub(1), col.saturating_sub(1)))
                        },
                        score: 0,
                        row,
                        col,
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        dp.iter_mut().enumerate().for_each(|(i, row)| {
            row[0].score = i as i32 * self.gap_score;
        });

        for i in 0..=a.len() {
            dp[0][i].score = i as i32 * self.gap_score;
        }

        for y in 1..dp.len() {
            for x in 1..dp[0].len() {
                self.fill_cell(&mut dp, y, x, a, b);
            }
        }

        let mut alignment: Vec<(Option<&T>, Option<&U>)> = vec![];
        let mut y = dp.len() - 1;
        let mut x = dp[y].len() - 1;
        let score = dp[y][x].score;
        while let Some(prev) = dp[y][x].prev {
            alignment.push((
                if dp[y][x].col - dp[prev.0][prev.1].col == 1 {
                    Some(&a[dp[y][x].col - 1])
                } else {
                    None
                },
                if dp[y][x].row - dp[prev.0][prev.1].row == 1 {
                    Some(&b[dp[y][x].row - 1])
                } else {
                    None
                },
            ));

            y = prev.0;
            x = prev.1;
        }

        alignment.reverse();
        Output { alignment, score }
    }

    fn fill_cell<T: PartialEq<U>, U: PartialEq<T>>(
        &self,
        dp: &mut [Vec<Cell>],
        y: usize,
        x: usize,
        a: &[T],
        b: &[U],
    ) {
        let row_space_score = dp[y - 1][x].score + self.gap_score;
        let col_space_score = dp[y][x - 1].score + self.gap_score;
        let mismatch_score = dp[y - 1][x - 1].score + self.replace(&b[y - 1], &a[x - 1]);

        if row_space_score >= col_space_score {
            if mismatch_score >= row_space_score {
                dp[y][x].score = mismatch_score;
                dp[y][x].prev = Some((y - 1, x - 1));
            } else {
                dp[y][x].score = row_space_score;
                dp[y][x].prev = Some((y - 1, x));
            }
        } else if mismatch_score >= col_space_score {
            dp[y][x].score = mismatch_score;
            dp[y][x].prev = Some((y - 1, x - 1));
        } else {
            dp[y][x].score = col_space_score;
            dp[y][x].prev = Some((y, x - 1));
        }
    }
}

#[derive(Default)]
struct Cell {
    prev: Option<(usize, usize)>,
    score: i32,
    row: usize,
    col: usize,
}

impl std::fmt::Debug for Cell {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.score.fmt(f)
    }
}

#[cfg(test)]
mod test {
    use super::Config;

    #[test]
    fn test() {
        let cases = [
            (
                "ACCCGGTCGTCAATTA",
                "ACCACCGGTTGGTCCAATAA",
                "A_C_CCGG_TCGT_CAATTA",
                "ACCACCGGTTGGTCCAATAA",
                8,
                Config::default(),
            ),
            (
                "ACCCGGTCGTCAATTA",
                "ACCACCGGTTGGTCCAATAA",
                "A_C_CCGG_T_CGT_CAAT_TA",
                "ACCACCGGTTG_GTCCAATA_A",
                6,
                Config::default().mismatch_score(-4),
            ),
            (
                "AGTACGCA",
                "TATGC",
                "AGTACGCA",
                "__TATGC_",
                1,
                Config {
                    match_score: 2,
                    mismatch_score: -1,
                    gap_score: -2,
                },
            ),
        ];

        for (a, b, expected_a, expected_b, expected_score, config) in cases {
            let a = a.chars().collect::<Vec<_>>();
            let b = b.chars().collect::<Vec<_>>();

            let output = config.compute(&a, &b);
            assert_eq!(output.score(), config.needleman_wunsch(&a, &b).score());

            let (aligned_a, aligned_b): (String, String) = output
                .alignment()
                .iter()
                .map(|(a, b)| (a.copied().unwrap_or('_'), b.copied().unwrap_or('_')))
                .unzip();

            assert_eq!(output.score(), expected_score);
            assert_eq!(
                [aligned_a, aligned_b],
                [expected_a.to_string(), expected_b.to_string()]
            );
        }
    }
}
