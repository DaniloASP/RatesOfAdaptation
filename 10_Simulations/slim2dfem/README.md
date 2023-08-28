Requires the bpp-core3 library from Bio++ V3 (earlier version will also work, but the makefile will need to be edited).
Adjust the CPPFLAGS and LDFLAGS if needed, based on where the library are installed.

```bash
make CPPFLAGS="-I$HOME/.local/include/" LDFLAGS="-L$HOME/.local/lib"
```

