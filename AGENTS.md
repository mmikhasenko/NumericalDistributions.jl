This repository uses Julia for testing. If `julia` is not available in the environment, install it with juliaup.

## Environment setup
1. Install Julia via juliaup:
   ```bash
   juliaup add 1
   juliaup default 1
   ```
2. Verify `julia --version` prints a valid version.

## Running tests
After code changes, run the test suite with:
```bash
julia --project -e 'using Pkg; Pkg.test()'
```
