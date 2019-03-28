# CallableLoci_processor

Combines CallableLoci bed files from control and experimental samples. Useful to determine the region in which mutations could be called. It generates bedfiles containing merged Callable regions for autosomes and for autosomes plus the X chromosome. It also calculates the total size of the Callable regions.

## USAGE:
```bash
python CallableLoci_processor.py dir_in dir_out exp_name control_sample
```
## Dependencies
- Python >= 3.6.5
- bedtools >= 2.27.1
- grep >= 3.1
- coreutils >= 8.30 (Standard unix tools)
