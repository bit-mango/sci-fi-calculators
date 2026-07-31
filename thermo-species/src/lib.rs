mod parsing;
mod types;
include!(concat!(env!("OUT_DIR"), "/generated_sample.rs"));

pub fn howdy() {
    // println!("Found a number in the generated file! {}", HELLO_WORLD);
}
