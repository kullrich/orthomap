Update/download NCBI database
=============================

You can download or update your local fo the NCBI's taxonomy database (~300MB)
using `orthomap`::

The following command downloads or updates your local copy of the
NCBI's taxonomy database (~300MB). The database is saved at
`~/.etetoolkit/taxa.sqlite`.

```python
>>> from orthomap import ncbitax
>>> ncbitax.update_ncbi()
```

`orthomap` uses the package 