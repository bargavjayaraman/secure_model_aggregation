Install obliv-c (https://github.com/samee/obliv-c) and its dependencies.

Install Cycle utility - https://github.com/samee/cmd

Compile using GCC wrapper of obliv-c
`/path/to/oblivcc secure_model_aggregation.c secure_model_aggregation.oc -I .`

Run: `cycle './a.out 1 | ./a.out 2'`

Aggregated model will be stored in Results folder along with timing and gate information.
