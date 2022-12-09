use crate::alphabet::{Alphabet, Scoring};
use crate::aminoacid::*;
use itertools::Itertools;
use std::fmt::Write;

/// An alignment of two reads.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// The score of this alignment
    pub score: isize,
    /// The path or steps taken for the alignment
    pub path: Vec<Piece>,
    /// The position in the first sequence where the alignment starts
    pub start_a: usize,
    /// The position in the second sequence where the alignment starts
    pub start_b: usize,
    /// The first sequence
    pub seq_a: Vec<AminoAcid>,
    /// The second sequence
    pub seq_b: Vec<AminoAcid>,
}

impl Alignment {
    fn short(&self) -> String {
        self.path.iter().map(Piece::short).join("")
    }

    fn aligned(&self) -> String {
        let blocks: Vec<char> = " ▁▂▃▄▅▆▇█".chars().collect();
        let blocks_neg: Vec<char> = " ▔▔▔▀▀▀▀█".chars().collect();
        let mut str_a = String::new();
        let mut str_b = String::new();
        let mut str_blocks = String::new();
        let mut str_blocks_neg = String::new();
        let mut loc_a = self.start_a;
        let mut loc_b = self.start_b;

        for piece in &self.path {
            let l = std::cmp::max(piece.step_b, piece.step_a);
            if piece.step_a == 0 {
                let _ = write!(str_a, "{:-<width$}", "", width = l as usize);
            } else {
                let _ = write!(
                    str_a,
                    "{:·<width$}",
                    self.seq_a[loc_a..loc_a + piece.step_a as usize]
                        .iter()
                        .map(std::string::ToString::to_string)
                        .join(""),
                    width = l as usize
                );
            }
            if piece.step_b == 0 {
                let _ = write!(str_b, "{:-<width$}", "", width = l as usize);
            } else {
                let _ = write!(
                    str_b,
                    "{:·<width$}",
                    self.seq_b[loc_b..loc_b + piece.step_b as usize]
                        .iter()
                        .map(std::string::ToString::to_string)
                        .join(""),
                    width = l as usize
                );
            }
            let _ = write!(
                str_blocks,
                "{}",
                str::repeat(
                    &if piece.local_score < 0 {
                        " ".to_string()
                    } else {
                        #[allow(clippy::cast_sign_loss)] // Checked above
                        blocks[piece.local_score as usize].to_string()
                    },
                    l as usize
                )
            );
            let _ = write!(
                str_blocks_neg,
                "{}",
                str::repeat(
                    &if piece.local_score > 0 {
                        " ".to_string()
                    } else {
                        #[allow(clippy::cast_sign_loss)] // Checked above
                        blocks_neg[-piece.local_score as usize].to_string()
                    },
                    l as usize
                )
            );

            loc_a += piece.step_a as usize;
            loc_b += piece.step_b as usize;
        }

        format!("{}\n{}\n{}\n{}", str_a, str_b, str_blocks, str_blocks_neg)
    }

    /// Generate a summary of this alignment for printing to the command line
    pub fn summary(&self) -> String {
        format!(
            "score: {}\npath: {}\nstart: ({}, {})\naligned:\n{}",
            self.score,
            self.short(),
            self.start_a,
            self.start_b,
            self.aligned()
        )
    }

    /// The total number of residues matched on the first sequence
    pub fn len_a(&self) -> usize {
        self.path.iter().map(|p| p.step_a as usize).sum()
    }

    /// The total number of residues matched on the second sequence
    pub fn len_b(&self) -> usize {
        self.path.iter().map(|p| p.step_b as usize).sum()
    }
}

/// A piece in an alignment, determining what step was taken in the alignment and how this impacted the score
#[derive(Clone, Default, Debug)]
pub struct Piece {
    /// The total score of the path up till now
    pub score: isize,
    /// The local contribution to the score of this piece
    pub local_score: i8,
    /// The number of steps on the first sequence
    pub step_a: u8,
    /// The number of steps on the second sequence
    pub step_b: u8,
}

impl Piece {
    /// Create a new alignment piece
    pub const fn new(score: isize, local_score: i8, step_a: u8, step_b: u8) -> Self {
        Self {
            score,
            local_score,
            step_a,
            step_b,
        }
    }
}

impl Piece {
    /// Display this piece very compactly
    pub fn short(&self) -> String {
        match (self.step_a, self.step_b) {
            (0, 1) => "I".to_string(),
            (1, 0) => "D".to_string(),
            (1, 1) => "M".to_string(),
            (a, b) => format!("S[{},{}]", b, a),
        }
    }
}

/// The type of alignment to perform
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Type {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    Global,
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    Local,
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    GlobalForB,
}

impl Type {
    const fn global(self) -> bool {
        !matches!(self, Self::Local)
    }
}

/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger then [`isize::MAX`].
#[allow(clippy::too_many_lines)]
pub fn align(seq_a: &[AminoAcid], seq_b: &[AminoAcid], alphabet: &Alphabet, ty: Type) -> Alignment {
    assert!(isize::try_from(seq_a.len()).is_ok());
    assert!(isize::try_from(seq_b.len()).is_ok());
    let mut matrix = vec![vec![Piece::default(); seq_b.len() + 1]; seq_a.len() + 1];
    let mut high = (0, 0, 0);

    if ty.global() {
        #[allow(clippy::cast_possible_wrap)]
        // b is always less than seq_b
        for index_b in 0..=seq_b.len() {
            matrix[0][index_b] = Piece::new(
                (index_b as isize) * Scoring::GapExtendPenalty as isize,
                Scoring::GapExtendPenalty as i8,
                0,
                if index_b == 0 { 0 } else { 1 },
            );
        }
    }
    if ty == Type::Global {
        #[allow(clippy::cast_possible_wrap)]
        // a is always less than seq_a
        for (index_a, row) in matrix.iter_mut().enumerate() {
            row[0] = Piece::new(
                (index_a as isize) * Scoring::GapExtendPenalty as isize,
                Scoring::GapExtendPenalty as i8,
                if index_a == 0 { 0 } else { 1 },
                0,
            );
        }
    }

    let mut values = Vec::with_capacity(Alphabet::STEPS * Alphabet::STEPS + 2);
    for index_a in 1..=seq_a.len() {
        for index_b in 1..=seq_b.len() {
            values.clear();
            for len_a in 0..=Alphabet::STEPS {
                for len_b in 0..=Alphabet::STEPS {
                    if len_a == 0 && len_b != 1
                        || len_a != 1 && len_b == 0
                        || len_a > index_a
                        || len_b > index_b
                    {
                        continue; // Do not allow double gaps (just makes no sense)
                    }
                    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                    // len_a and b are always less then Alphabet::STEPS
                    let score = if len_a == 0 || len_b == 0 {
                        Scoring::GapExtendPenalty as i8 // Defined to always be one gap
                    } else {
                        alphabet[(
                            &seq_a[index_a - len_a..index_a],
                            &seq_b[index_b - len_b..index_b],
                        )]
                    };
                    if score == 0 {
                        continue;
                    }
                    values.push(Piece::new(
                        matrix[index_a - len_a][index_b - len_b].score + score as isize,
                        score,
                        len_a as u8,
                        len_b as u8,
                    ));
                }
            }
            let value = values
                .iter()
                .max_by(|x, y| x.score.cmp(&y.score))
                .cloned()
                .unwrap_or_default();
            if value.score >= high.0 {
                high = (value.score, index_a, index_b);
            }
            matrix[index_a][index_b] = value;
        }
    }

    // loop back
    if ty == Type::Global {
        high = (
            matrix[seq_a.len()][seq_b.len()].score,
            seq_a.len(),
            seq_b.len(),
        );
    } else if ty == Type::GlobalForB {
        let value = (0..=seq_a.len())
            .map(|v| (v, matrix[v][seq_b.len()].score))
            .max_by(|a, b| a.1.cmp(&b.1))
            .unwrap_or_default();
        high = (value.1, value.0, seq_b.len());
    }
    let mut path = Vec::new();
    let high_score = high.0;
    //dbg!(&highest_score);
    //dbg!(&matrix);
    while !(high.1 == 0 && high.2 == 0) {
        let value = matrix[high.1][high.2].clone();
        if value.step_a == 0 && value.step_b == 0 {
            break;
        }
        high = (
            0,
            high.1 - value.step_a as usize,
            high.2 - value.step_b as usize,
        );
        path.push(value);
    }
    //dbg!(&path);
    Alignment {
        score: high_score,
        path: path.into_iter().rev().collect(),
        start_a: high.1,
        start_b: high.2,
        seq_a: seq_a.to_owned(),
        seq_b: seq_b.to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::{align, Type};
    use crate::alphabet::Alphabet;
    use crate::aminoacid::AminoAcid::*;

    #[test]
    fn equal() {
        let alphabet = Alphabet::default();
        let a = vec![A, C, C, G, W];
        let b = vec![A, C, C, G, W];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        assert_eq!(40, result.score);
        assert_eq!("MMMMM", &result.short());
    }

    #[test]
    fn insertion() {
        let alphabet = Alphabet::default();
        let a = vec![A, C, G, W];
        let b = vec![A, C, F, G, W];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        assert_eq!(27, result.score);
        assert_eq!("MMIMM", &result.short());
    }

    #[test]
    fn deletion() {
        let alphabet = Alphabet::default();
        let a = vec![A, C, F, G, W];
        let b = vec![A, C, G, W];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        assert_eq!(27, result.score);
        assert_eq!("MMDMM", &result.short());
    }

    #[test]
    fn iso_mass() {
        let alphabet = Alphabet::default();
        let a = vec![A, F, G, G, W];
        let b = vec![A, F, N, W];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        dbg!(result.short());
        assert_eq!(29, result.score);
        assert_eq!("MMS[1,2]M", &result.short());
    }

    #[test]
    fn switched() {
        let alphabet = Alphabet::default();
        let a = vec![A, F, G, G, W];
        let b = vec![A, G, F, G, W];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        dbg!(result.short());
        assert_eq!(28, result.score);
        assert_eq!("MS[2,2]MM", &result.short());
    }

    #[test]
    fn local() {
        let alphabet = Alphabet::default();
        let a = vec![A, F, G, G, E, W];
        let b = vec![F, G, G, D];
        let result = align(&a, &b, &alphabet, Type::Local);
        dbg!(&result);
        dbg!(result.short());
        assert_eq!(24, result.score);
        assert_eq!("MMM", &result.short());
    }

    #[test]
    fn global() {
        let alphabet = Alphabet::default();
        let a = vec![A, F, G, G, E, W];
        let b = vec![F, G, G, D];
        let result = align(&a, &b, &alphabet, Type::Global);
        dbg!(&result);
        println!("{}", result.summary());
        assert_eq!(13, result.score);
        assert_eq!("DMMMDM", &result.short());
        assert_eq!(0, result.start_a, "A global alignment should start at 0");
    }

    #[test]
    fn global_for_b() {
        let alphabet = Alphabet::default();
        let a = vec![A, F, G, G, E, W];
        let b = vec![F, G, G, D];
        let result = align(&a, &b, &alphabet, Type::GlobalForB);
        dbg!(&result);
        dbg!(result.short());
        assert_eq!(23, result.score);
        assert_eq!("MMMM", &result.short());
        assert_eq!(0, result.start_b, "A global alignment should start at 0");
    }
}
