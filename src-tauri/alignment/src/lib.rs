//! An algorithm based on Needleman Wunsch/Smith Waterman but extended to allow for mass based alignment.
//! The mass based part gives the option to match two sets of aminoacids with different sizes.
//! For example the set {Q} matches {AG} because these have the same mass and so are commonly misclassified
//! is de novo sequencing for peptides. Besides iso mass definitions it also handles swaps with finesse,
//! meaning that {AG} matches {GA} with a well defined score to allow for these mistakes to be fixed. The
//! last important addition is the handling of post translational modifications meaning that {Q} matches {E}
//! but not the other way around to allow for deamidation of the sample in reference to the template.
//!
//! ```rust
//! use mass_alignment::*;
//! use mass_alignment::AminoAcid::*;
//!
//! let alphabet = Alphabet::default();
//! let template = &[A,G,Q,S,T,Q];
//! let query = &[Q,E,S,W];
//! let result = align(template, query, &alphabet, Type::GlobalForB);
//! println!("{}", result.summary());
//! assert_eq!(15, result.score)
//! ```

#![allow(dead_code)]
#![warn(clippy::pedantic, clippy::nursery, clippy::all, missing_docs)]
#![allow(
    clippy::enum_glob_use,
    clippy::wildcard_imports,
    clippy::must_use_candidate
)]
/// The module containing all alignment handling
mod alignment;
/// The module containing all alphabet handling
mod alphabet;
/// The module containing the definition for aminoacids
mod aminoacid;
/// The module containing the definition for templates
pub mod template;

pub use crate::alignment::align;
pub use crate::alignment::*;
pub use crate::alphabet::Alphabet;
pub use crate::aminoacid::sequence_from_string;
pub use crate::aminoacid::sequence_to_string;
pub use crate::aminoacid::AminoAcid;
