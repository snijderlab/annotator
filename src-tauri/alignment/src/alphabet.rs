use crate::aminoacid::AminoAcid;
use crate::aminoacid::AminoAcid::*;
use itertools::Itertools;

/// An alphabet to determine the score of two amino acid sets
pub struct Alphabet {
    array: Vec<Vec<i8>>,
}

impl std::ops::Index<(&[AminoAcid], &[AminoAcid])> for Alphabet {
    type Output = i8;
    fn index(&self, index: (&[AminoAcid], &[AminoAcid])) -> &Self::Output {
        &self.array[get_index(index.0)][get_index(index.1)]
    }
}

fn get_index_ref(set: &[&AminoAcid]) -> usize {
    set.iter()
        .fold(0, |acc, item| acc * AminoAcid::MAX + **item as usize)
}

fn get_index(set: &[AminoAcid]) -> usize {
    set.iter()
        .fold(0, |acc, item| acc * AminoAcid::MAX + *item as usize)
}

impl Alphabet {
    /// The number of steps to trace back, if updated a lot of other code has to be updated as well
    pub const STEPS: usize = 3;
}

#[repr(i8)]
#[derive(Clone, Default, Debug)]
pub enum Scoring {
    /// The score for identity, should be the highest score of the bunch
    Identity = 8,
    /// The score for a mismatch
    #[default]
    Mismatch = -1,
    /// The score for an iso mass set, eg Q<>AG
    IsoMass = 5,
    /// The score for a modification
    Modification = 3,
    /// The score for a switched set, defined as this value times the size of the set (eg AG scores 4 with GA)
    Switched = 2,
    /// The score for scoring a gap, should be less than `MISMATCH`
    GapStartPenalty = -5,
    GapExtendPenalty = -3,
}

#[allow(clippy::too_many_lines)]
impl Default for Alphabet {
    fn default() -> Self {
        macro_rules! sets {
            ($($($($id:ident),+);+)|+) => {
                vec![
                    $(vec![
                        $(vec![$($id),+],)+
                      ],)+
                ]
            };
        }

        #[allow(clippy::cast_possible_truncation)]
        // STEPS is always within bounds for u32
        let mut alphabet = Self {
            array: vec![
                vec![0; (AminoAcid::MAX + 1).pow(Self::STEPS as u32)];
                (AminoAcid::MAX + 1).pow(Self::STEPS as u32)
            ],
        };

        for x in 0..=AminoAcid::MAX {
            for y in 0..=AminoAcid::MAX {
                alphabet.array[x][y] = if x == y {
                    Scoring::Identity as i8
                } else {
                    Scoring::Mismatch as i8
                };
            }
        }
        let iso_mass = sets!(
            I;     L|
            N;     G,G|
            Q;     A,G|
            A,V;   G,L;   G,I|
            A,N;   Q,G;   A,G,G|
            L,S;   I,S;   T,V|
            A,M;   C,V|
            N,V;   A,A,A; G,G,V|
            N,T;   Q,S;   A,G,S; G,G,T|
            L,N;   I,N;   Q,V;   A,G,V; G,G,L; G,G,I|
            D,L;   D,I;   E,V|
            Q,T;   A,A,S; A,G,T|
            A,Y;   F,S|
            L,Q;   I,Q;   A,A,V; A,G,L; A,G,I|
            N,Q;   A,N,G; Q,G,G|
            K,N;   G,G,K|
            E,N;   D,Q;   A,D,G; E,G,G|
            D,K;   A,A,T; G,S,V|
            M,N;   A,A,C; G,G,M|
            A,S;   G,T|
            A,A,L; A,A,I; G,V,V|
            Q,Q;   A,A,N; A,Q,G|
            E,Q;   A,A,D; A,E,G|
            E,K;   A,S,V; G,L,S; G,I,S; G,T,V|
            M,Q;   A,G,M; C,G,V|
            A,A,Q; N,G,V
        );

        for set in iso_mass {
            for set in set.iter().permutations(2) {
                let a = set[0];
                let b = set[1];
                for seq_a in a.iter().permutations(a.len()) {
                    for seq_b in b.iter().permutations(b.len()) {
                        alphabet.array[get_index_ref(&seq_a)][get_index_ref(&seq_b)] =
                            Scoring::IsoMass as i8;
                    }
                }
            }
        }

        let modifications = sets!(
            //N;D| // Amidation only at N term
            Q;E| // Deamidation
            D;N| // Deamidation
            C;T| // Disulfide bond
            T;D| // Methylation
            S;T| // Methylation
            D;E| // Methylation
            R;A,V;G,L| // Methylation
            Q;A,A // Methylation
        );

        for set in modifications {
            let a = &set[0];
            for seq_b in set.iter().skip(1) {
                alphabet.array[get_index(a)][get_index(seq_b.as_slice())] =
                    Scoring::Modification as i8;
            }
        }

        let amino_acids = (1..=AminoAcid::MAX)
            .map(|a| AminoAcid::try_from(a).unwrap())
            .collect_vec();
        for size in 2..=Self::STEPS {
            for set in amino_acids
                .iter()
                .combinations_with_replacement(size)
                .flat_map(|v| v.into_iter().permutations(size))
            {
                if set.iter().all(|v| *v == set[0]) {
                    continue; // Do not add [A, A] or [A, A, A] etc as SWITCHED
                }
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                // set.len() is at max equal to Self::STEPS
                for switched in set.clone().into_iter().permutations(size) {
                    alphabet.array[get_index_ref(&set)][get_index_ref(&switched)] =
                        Scoring::Switched as i8 * set.len() as i8;
                }
            }
        }

        alphabet
    }
}

#[cfg(test)]
mod tests {
    use super::{Alphabet, Scoring};
    use crate::aminoacid::AminoAcid::*;

    #[test]
    fn identity() {
        let alphabet = Alphabet::default();
        assert_eq!(
            Scoring::Identity as i8,
            alphabet[([A].as_slice(), [A].as_slice())]
        );
        assert_eq!(0, alphabet[([A, A].as_slice(), [A, A].as_slice())]);
        assert_eq!(0, alphabet[([A, A, A].as_slice(), [A, A, A].as_slice())]);
    }

    #[test]
    fn similarity() {
        let alphabet = Alphabet::default();
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([I].as_slice(), [L].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([N].as_slice(), [G, G].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([G, G].as_slice(), [N].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([A, S].as_slice(), [G, T].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([S, A].as_slice(), [G, T].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([S, A].as_slice(), [T, G].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([A, S].as_slice(), [T, G].as_slice())]
        );
        assert_eq!(
            Scoring::IsoMass as i8,
            alphabet[([L, Q].as_slice(), [A, V, A].as_slice())]
        );
    }

    #[test]
    fn inequality() {
        let alphabet = Alphabet::default();
        assert_eq!(
            Scoring::Mismatch as i8,
            alphabet[([I].as_slice(), [Q].as_slice())]
        );
        assert_eq!(0, alphabet[([Q].as_slice(), [G, G].as_slice())]);
        assert_eq!(0, alphabet[([A, E].as_slice(), [G, T].as_slice())]);
        assert_eq!(0, alphabet[([E, Q].as_slice(), [A, V, A].as_slice())]);
    }

    #[test]
    fn switched() {
        let alphabet = Alphabet::default();
        assert_eq!(
            Scoring::Switched as i8 * 2,
            alphabet[([E, Q].as_slice(), [Q, E].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([D, A, C].as_slice(), [A, C, D].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([C, D, A].as_slice(), [A, C, D].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([A, C, D].as_slice(), [D, A, C].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([C, D, A].as_slice(), [D, A, C].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([A, C, D].as_slice(), [C, D, A].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([D, A, C].as_slice(), [C, D, A].as_slice())]
        );
        assert_eq!(
            Scoring::Switched as i8 * 3,
            alphabet[([V, A, A].as_slice(), [A, V, A].as_slice())]
        );
    }

    #[test]
    fn modification() {
        let alphabet = Alphabet::default();
        assert_eq!(
            Scoring::Modification as i8,
            alphabet[([D].as_slice(), [N].as_slice())]
        );
        assert_eq!(
            Scoring::Mismatch as i8,
            alphabet[([N].as_slice(), [D].as_slice())]
        );
    }
}
