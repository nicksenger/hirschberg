# Hirschberg

Hirschberg provides a generic implementation of the Hirschberg algorithm for finding the optimal global alignment of two sequences. 

Example usage:

```rust
let a = "ACCCGGTCGTCAATTA".chars().collect::<Vec<_>>();
let b = "ACCACCGGTTGGTCCAATAA".chars().collect::<Vec<_>>();

let output = hirschberg::Config::default().compute(&a, &b);

let (aligned_a, aligned_b): (String, String) = output
    .alignment()
    .iter()
    .map(|[a, b]| (a.copied().unwrap_or('_'), b.copied().unwrap_or('_')))
    .unzip();

assert_eq!(output.score(), 8);
assert_eq!(
    [aligned_a, aligned_b],
    [
        "A_C_CCGG_TCGT_CAATTA".to_string(),
        "ACCACCGGTTGGTCCAATAA".to_string()
    ]
);
```

The match, mismatch and gap contribution to the distance may be configured:

```rust
let a = "ACCCGGTCGTCAATTA".chars().collect::<Vec<_>>();
let b = "ACCACCGGTTGGTCCAATAA".chars().collect::<Vec<_>>();

let output = hirschberg::Config::default().mismatch_score(-4).compute(&a, &b);

let (aligned_a, aligned_b): (String, String) = output
    .alignment()
    .iter()
    .map(|[a, b]| (a.copied().unwrap_or('_'), b.copied().unwrap_or('_')))
    .unzip();

assert_eq!(output.score(), 6);
assert_eq!(
    [aligned_a, aligned_b],
    [
        "A_C_CCGG_T_CGT_CAAT_TA".to_string(),
        "ACCACCGGTTG_GTCCAATA_A".to_string()
    ]
);
```
