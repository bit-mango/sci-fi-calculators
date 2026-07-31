Phase 1 — Scaffold the crate
1. Decide a name now even though you won't publish yet (crates.io names are first-come; check availability so you're not stuck retitling later). Something like thermo-species or nasa-thermo.
2. Create it as a new workspace member (cargo new <name> --lib), add it to the root Cargo.toml's members list alongside fusion-engine and mission-design.
3. Move thermo.inp into the new crate's own directory (e.g. a data/ subfolder inside it) — not left in fusion-engine.
4. Fill in real [package] fields in its Cargo.toml now (license, description, a starting version = "0.1.0") even though nothing depends on it yet.
5. Confirm it builds empty (cargo build -p <name>) before touching build.rs — establishes your baseline.

Phase 2 — Get oriented on build.rs mechanics (no parsing yet)
6. Add a build.rs at the crate root (same level as Cargo.toml, not in src/). Cargo auto-detects and runs it before compiling src/.
7. Learn the two things build.rs communicates through: (a) it can write files to the path in the OUT_DIR env var — that's where your generated Rust code goes; (b) it talks back to cargo via cargo: prefixed lines printed to stdout — you'll want cargo:rerun-if-changed=<path> for your input file and the env var, so cargo only re-runs the build script when the data actually changes.
8. Do a trivial dry run first: have build.rs write a tiny hardcoded .rs file into OUT_DIR (e.g. one const), and include!() it from src/lib.rs. Get this loop working end to end before any real parsing — it's the "hello world" of codegen and isolates build.rs mechanics from parsing bugs.

Phase 3 — Parse one record by hand first
9. Before writing the general parser, manually work out the field layout for a single species block using the e- and Ag entries you already looked at — write down (on paper or a comment) exactly which whitespace-delimited token on which line maps to which field. Don't start coding the parser until you can say precisely "token 3 of line 2 is the phase flag" etc.
10. Write the parser to walk the file line-by-line, using a small state machine: "expecting header line" → "expecting metadata line" → "expecting interval line 1 of 3" → "expecting coefficient line 1 of 2" → "expecting coefficient line 2 of 2" → loop back based on how many intervals the metadata line declared.
11. Handle the Fortran D exponent notation (e.g. -7.453750000D+02) — you'll need to swap D/d for E before f64::parse, since Rust's parser doesn't recognize D.
12. Test the parser against just the first 2-3 species in the file before running it on all ~2000 — print/assert the parsed struct matches what you hand-decoded in step 9.

Phase 4 — Design the generated output
13. Decide the Species struct shape (fields: symbol, MW, phase, ΔH_f, temperature ranges each with their 9 coefficients, element composition) — this is worth sketching on paper before generating any code, since the codegen just needs to print Rust syntax for whatever struct you land on.
14. Write the codegen step: for each parsed species, emit one pub const <NAME>: Species = Species { ... }; inside an impl Species block, plus a from_str match arm. Sanitize names into valid Rust identifiers (species names contain characters like +, -, parentheses that aren't valid in identifiers — decide your mangling scheme now, e.g. Ag+ → AgPlus).
15. Run the full file through and just check it compiles — don't worry about correctness of every species yet, just that codegen produces valid Rust for all ~2000 entries.

Phase 5 — Env var override for custom files
16. Add logic in build.rs to check for your chosen env var (e.g. THERMO_INP_PATH) and fall back to the bundled data/thermo.inp if unset. Add cargo:rerun-if-env-changed=THERMO_INP_PATH so cargo knows to reparse when it's set/changed.

Phase 6 — Feature-gate by phase
17. Since the phase flag is already in the parsed data, add gas / condensed (and default = both, or pick a sensible default) features to Cargo.toml, and wrap each generated const in the matching #[cfg(feature = "...")] at codegen time based on that species' phase.
18. Test building with --no-default-features --features gas and confirm only gas-phase consts exist (a condensed-only reference should fail to compile).

Phase 7 — Wire into fusion-engine
19. Add the new crate as a path dependency in fusion-engine's Cargo.toml.
20. Migrate fusion-engine's existing hand-written Species enum usages over to the new crate incrementally — don't do a big-bang swap; pick one call site, confirm it compiles and behaves the same, then move to the next.

Once this is working end to end, that's a natural point to circle back on the species!() macro and the fine-grained allowlist, per what you already deferred. Ask me anything as you go — happy to explain build.rs quirks (like why OUT_DIR moves between debug/release, or println! vs cargo: output) or sanity-check your field-mapping from step 9 before you commit to a struct shape.
