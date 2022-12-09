use crate::alignment::Alignment;
use crate::alphabet::Scoring;
use crate::aminoacid::{self, AminoAcid};
use crate::{align, Alphabet};
use itertools::Itertools;
use std::fmt::Write;

/// A template that is matched with many reads
pub struct Template {
    /// The sequence of this template
    pub sequence: Vec<AminoAcid>,
    /// The reads matched to this template
    pub reads: Vec<Alignment>,
}

impl Template {
    /// Create a new template by matching the given reads to the given template sequence
    pub fn new(sequence: Vec<AminoAcid>, reads: Vec<&[AminoAcid]>, alphabet: &Alphabet) -> Self {
        Self {
            reads: reads
                .into_iter()
                .map(|v| align(&sequence, v, alphabet, crate::alignment::Type::GlobalForB))
                .collect(),
            sequence,
        }
    }

    /// Generate HTML for a reads alignment, all styling is missing and it is only a small part of a document
    pub fn generate_html(&self) -> String {
        let mut insertions = vec![0; self.sequence.len()];
        for read in &self.reads {
            let mut loc_a = read.start_a;
            let mut insertion = 0;

            for piece in &read.path {
                if piece.step_a == 0 && piece.step_b == 1 {
                    insertion += 1;
                } else if insertion != 0 {
                    insertions[loc_a] = std::cmp::max(insertions[loc_a], insertion);
                    insertion = 0;
                } else {
                    insertion = 0;
                }
                loc_a += piece.step_a as usize;
            }
        }

        let mut output = format!("<div class='alignment-body' style='grid-template-columns: repeat({}, 1ch);'><div class='numbering'>1....</div><div class='template'>", insertions.iter().sum::<usize>() + self.sequence.len());
        for (ins, seq) in insertions.iter().zip(&self.sequence) {
            let _ = write!(output, "{:-<1$}{2}", "", ins, seq);
        }
        let _ = write!(output, "</div>");

        for read in &self.reads {
            let _ = write!(
                output,
                "<div style='grid-column: {} / {}'>",
                insertions[0..read.start_a].iter().sum::<usize>() + read.start_a + 1,
                insertions[0..read.start_a + read.len_a()]
                    .iter()
                    .sum::<usize>()
                    + read.start_a
                    + read.len_a()
                    + 3
            );
            let mut loc_a = read.start_a;
            let mut loc_b = read.start_b;
            let mut insertion = 0;
            for piece in &read.path {
                if piece.step_a == 0 && piece.step_b == 1 {
                    insertion += 1;
                } else {
                    let _ = write!(output, "{:-<1$}", "", insertions[loc_a] - insertion);
                    insertion = 0;
                }
                let _ = write!(
                    output,
                    "{}",
                    match (piece.step_a, piece.step_b) {
                        (0 | 1, 1) => read.seq_b[loc_b].to_string(),
                        (1, 0) => "-".to_string(),
                        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                        // a is defined to be in range 0..=Alphabet::STEPS
                        (a, b) => {
                            let inner = if a == b {
                                // As a equals b it is a swap or iso length iso mass sets, add in the missing insertions (if any)
                                // Because a equals b the length of the sequence patch and insertion patch is always equal.
                                // This means that the resulting insertions makes the text nicely aligned.
                                read.seq_b[loc_b..loc_b + b as usize]
                                    .iter()
                                    .zip(&insertions[loc_a..loc_a + a as usize])
                                    .map(|(sb, sa)| format!("{:->1$}", sb.to_string(), sa + 1))
                                    .join("")
                            } else {
                                aminoacid::sequence_to_string(
                                    &read.seq_b[loc_b..loc_b + b as usize],
                                )
                            };
                            format!(
                                "<span class='s{}' style='--i:{};--w:{};'>{}</span>",
                                if a == b && piece.local_score == Scoring::Switched as i8 * a as i8
                                {
                                    " swap"
                                } else {
                                    ""
                                },
                                inner.len(),
                                insertions[loc_a..loc_a + a as usize].iter().sum::<usize>()
                                    + a as usize,
                                inner
                            )
                        }
                    }
                );
                loc_a += piece.step_a as usize;
                loc_b += piece.step_b as usize;
            }
            let _ = write!(output, "</div>");
        }
        let _ = write!(output, "</div>");
        output
    }
}
