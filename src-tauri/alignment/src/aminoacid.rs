use itertools::Itertools;
use std::fmt::Display;

/// All aminoacids
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Eq, Ord)]
pub enum AminoAcid {
    /// Alanine
    A = 1,
    /// Arginine
    R,
    /// Asparagine
    N,
    /// Aspartic acid
    D,
    /// Cysteine
    C,
    /// Glutamine
    Q,
    /// Glutamic acid
    E,
    /// Glycine
    G,
    /// Histidine
    H,
    /// Isoleucine
    I,
    /// Leucine
    L,
    /// Lysine
    K,
    /// Methionine
    M,
    /// Phenylalanine
    F,
    /// Proline
    P,
    /// Serine
    S,
    /// Threonine
    T,
    /// Tryptophan
    W,
    /// Tyrosine
    Y,
    /// Valine
    V,
    /// Weird
    B,
    /// Also weird
    Z,
    /// Single gap
    X,
    /// Longer gap
    Gap,
}

impl AminoAcid {
    /// The total number of normal amino acids (disregards Gap)
    pub const MAX: usize = 23;
}

impl Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::A => "A",
            Self::R => "R",
            Self::N => "N",
            Self::D => "D",
            Self::C => "C",
            Self::Q => "Q",
            Self::E => "E",
            Self::G => "G",
            Self::H => "H",
            Self::I => "I",
            Self::L => "L",
            Self::K => "K",
            Self::M => "M",
            Self::F => "F",
            Self::P => "P",
            Self::S => "S",
            Self::T => "T",
            Self::W => "W",
            Self::Y => "Y",
            Self::V => "V",
            Self::B => "B",
            Self::Z => "Z",
            Self::X => "X",
            Self::Gap => "*",
        })
    }
}

impl TryFrom<usize> for AminoAcid {
    type Error = ();
    fn try_from(num: usize) -> Result<Self, Self::Error> {
        match num {
            1 => Ok(Self::A),
            2 => Ok(Self::R),
            3 => Ok(Self::N),
            4 => Ok(Self::D),
            5 => Ok(Self::C),
            6 => Ok(Self::Q),
            7 => Ok(Self::E),
            8 => Ok(Self::G),
            9 => Ok(Self::H),
            10 => Ok(Self::I),
            11 => Ok(Self::L),
            12 => Ok(Self::K),
            13 => Ok(Self::M),
            14 => Ok(Self::F),
            15 => Ok(Self::P),
            16 => Ok(Self::S),
            17 => Ok(Self::T),
            18 => Ok(Self::W),
            19 => Ok(Self::Y),
            20 => Ok(Self::V),
            21 => Ok(Self::B),
            22 => Ok(Self::Z),
            23 => Ok(Self::X),
            24 => Ok(Self::Gap),
            _ => Err(()),
        }
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = ();
    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'A' => Ok(Self::A),
            'R' => Ok(Self::R),
            'N' => Ok(Self::N),
            'D' => Ok(Self::D),
            'C' => Ok(Self::C),
            'Q' => Ok(Self::Q),
            'E' => Ok(Self::E),
            'G' => Ok(Self::G),
            'H' => Ok(Self::H),
            'I' => Ok(Self::I),
            'L' => Ok(Self::L),
            'K' => Ok(Self::K),
            'M' => Ok(Self::M),
            'F' => Ok(Self::F),
            'P' => Ok(Self::P),
            'S' => Ok(Self::S),
            'T' => Ok(Self::T),
            'W' => Ok(Self::W),
            'Y' => Ok(Self::Y),
            'V' => Ok(Self::V),
            'B' => Ok(Self::B),
            'Z' => Ok(Self::Z),
            'X' => Ok(Self::X),
            '*' => Ok(Self::Gap),
            _ => Err(()),
        }
    }
}

/// Create an aminoacid sequence from a string, just ignores any non aminoacids characters
pub fn sequence_from_string(value: &str) -> Vec<AminoAcid> {
    value
        .chars()
        .filter_map(|v| AminoAcid::try_from(v).ok())
        .collect()
}

/// Generate a string from a sequence of aminoacids
pub fn sequence_to_string(value: &[AminoAcid]) -> String {
    value.iter().map(std::string::ToString::to_string).join("")
}
