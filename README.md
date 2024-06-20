* Linear Stochastic estimate code that was developed my PhD in 2013.

This code computes the linear estimate coefficients ($L_{ij}$) to compute conditional eddies. This is computed using equation 2.2 in document (thesis) https://doi.org/10.4233/uuid:e9962229-d045-4614-876e-de3e7e5f188f

* Usage
** Requirements
   - Intel compiler
   - MKL libraries for FFTWs and LAPACK
   - GCC can also be used with LAPACK, FFTW

```bash
	$make      # create excutable
	$./corr    # compute LSE correlations
```

** Settings 
   In main file one can edit
   - first file to read
   - last file to read
   - number of files
   - event plane location needs to be set

   In mod_geometry :
   - fix parameters connected with size of domain and all
   - path to database
   - path for output



