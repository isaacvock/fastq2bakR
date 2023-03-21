extern crate rust_htslib;

use rust_htslib::bam::{self, Read};

fn count_tc_mismatches(bam_path: &str) -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    let mut bam = bam::Reader::from_path(bam_path)?;
    let header = bam.header().clone();
    let mut tc_mismatches = vec![];

    for r in bam.records() {
        let record = r?;
        let seq = record.seq();
        let qual = record.qual();
        let mut read_tc_mismatches = 0;

        for i in 0..seq.len() {
            if seq[i] == b'T' && record.is_reverse() {
                let ref_base = header.reference_id2name(record.tid()).unwrap()[record.pos() as usize] as char;

                if ref_base == 'C' && qual[i] >= 30 {
                    read_tc_mismatches += 1;
                }
            }
        }

        tc_mismatches.push(read_tc_mismatches);
    }

    Ok(tc_mismatches)
}
