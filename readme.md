# ERHAMB

ERHAMB is a modified version of Zbigniew Kisiel's ERHAMZ Version V16g-R3 which in turn is a modified version of Peter Groner's ERHAM Version V16g-R3.

Both versions are provided on the [PROSPE webpage](http://www.ifpan.edu.pl/~kisiel). For citing ERHAM please use the  references stated on the [PROSPE webpage](http://www.ifpan.edu.pl/~kisiel): P. Groner, *J. Chem. Phys.* **107**, 4483-4498 (1997) and P. Groner, *J. Mol.Spectrosc.* **278**, 52-67 (2012).

The here made changes increase the array sizes to allow for larger datasets and make it compilable with the gfortran compiler.

Compile with:

```bash
gfortran -O2 -w -std=legacy -ffixed-line-length-82 ./erhamb.for -o erhamb
```


