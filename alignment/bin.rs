#![allow(dead_code)]
#![warn(clippy::pedantic, clippy::nursery, clippy::all)]
#![allow(clippy::enum_glob_use, clippy::wildcard_imports)]
use mass_alignment::template::Template;
use mass_alignment::*;

fn main() {
    let alphabet = Alphabet::default();
    //let template = aminoacid::sequence_from_string("XXXXXXXXXXXXXXXXXXXXYFDYWGQGTLVTVSS");
    let template = mass_alignment::sequence_from_string("EVQLVESGGGLVQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARXXXXXXXXXXXXXXXXXXXX");
    let reads: Vec<Vec<AminoAcid>> = [
        //"SRWGGDGFYAMDYWGQGTLVTV",
        //"DWNGFYAMDYWGQGTLVTVSS",
        //"RWGGDGFYAMDYWGQGTLVTV",
        //"HVPHGDGFYAMDYWGQGTLVT",
        //"WRGGDGFYAMDYWGQGTLVT",
        //"SRWGGDGFYAMDYWGQGTLV",
        //"RWGGDGFYAMDYWGQGTLVT",
        //"WRNDGFYAMDYWGQGTLVT",
        //"RWGGDGFYAMDYWGQGTLV",
        //"MARNDGFYAMDYWGQGTLV",
        //"RWNDGFYAMDYWGQGTLV",
        //"SRWGGNGFYWDYWGQGT",
        //"RWNDGFYWDYWGQGT",
        //"DYWGQGTLVVTSS",
        //"DYWGQGTLVTVSS",
        //"DYWGQGTLVTV",
        //"DYWGQGTLVT",
        //"WGQGTLVT",
        "DLQLVESGGGLVGAKSPPGTLSAAASGFNL",
        "DLQLVESGGGLVGAKSPPGTLSAAASGFNL",
        "EVQLVESGGGLVQPGGSLSGAKYHSGFNL",
        "EVVQLVESGGGLVQPGGSLGVLSCAASGF",
        "DLQLVESGGGLVQPGGSLGVLSCAASGF",
        "DLQLVESGGGLVQPGTPLYWNAASGFNL",
        "DLQLVESGGGLVQPGGSLRLSCAASGF",
        "QVQLVESGGGLVQPGGSLRLSCAASGF",
        "EVQLVESGGGLPVQGGSLRLSCAADGF",
        "EVQLVESGGGLVQPGGSLRLSCAASGF",
        "EVQLVSGEGGLVQPGGSLRLSCAASGF",
        "QVELVESGGGLVQPGGSLRLSCAASGF",
        "TLSADTSKNTAYLQMNSLRAEDTAVY",
        "RFTLSADTSKNTAYLQMNSLRAEDTA",
        "QLVESGGGLVQPGGSLTHVAGAGHSGF",
        "SADTSKNTAYLQMNSLRAEDTAVYY",
        "LMLTDGYTRYADSVKGRFTLSADTS",
        "QLVESGGGLVQPGGSLRLSCAASGF",
        "QLVESGGGLVQPGGSLRLSCQTGF",
        "LVESGGGLVQPNSLRLSCAASGF",
    ]
    .into_iter()
    .map(mass_alignment::sequence_from_string)
    .collect();

    let template = Template::new(
        template,
        reads.iter().map(std::vec::Vec::as_slice).collect(),
        &alphabet,
    );
    let content = format!(
        "<!DOCTYPE html>
<html>
<head>
    <link rel='stylesheet' href='C:/Users/5803969/src/stitch/assets/styles.css'>
    <style>
        .s:not(.swap) {{
            width: calc(var(--w) * 1.001ch);
            transform:
                scaleX(calc(var(--w) / var(--i)));
            display: inline-block;
            color: var(--color-halfway);
            transform-origin: center left;
        }}

        .s.swap,
        .s.swap {{
            text-decoration: underline var(--color-halfway);
        }}

        .alignment .alignment-body .template,
        .alignment .alignment-body .numbering {{
            grid-column: 1/ -1;
        }}

        .alignment .alignment-wrapper {{
            margin-left: 0;
        }}
    </style>
</head>

<body>
    <div class='alignment'>
        <div class='alignment-wrapper'>{}
        </div>
    </div>
</body>    
</html>",
        template.generate_html()
    );
    std::fs::write("test.html", content).unwrap();
}
