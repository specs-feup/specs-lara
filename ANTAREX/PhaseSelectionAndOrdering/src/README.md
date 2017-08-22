# Compiler Phase Selection and Phase Ordering Using GCC and LLVM

## Authors
Ricardo Nobre   - ricardo.nobre@fe.up.pt
https://web.fe.up.pt/~pro11012/index.php?lang=en

João Bispo      - jbispo@fe.up.pt
http://jbispo.eu/

Tiago Carvalho	- tiago.carvalho@fe.up.pt


## License
This software is published under the GNU Lesser General Public License (LGPL) version 3.

## Requirements

### Clava
```
source: https://github.com/specs-feup/clava/tree/master/ClavaWeaver
binary: http://specs.fe.up.pt/tools/clava.jar
```

## Examples

### Using LLVM

#### Minimum options (using defaults)
```
clava ../src/LaradLauncher.lara -av "('../src', 'kernel_nussionov.c', 'llvm')" -nw -nci -b 2
```
#### All options (with named parameters)
```
clava ../src/LaradLauncher.lara -av "{laradFoldername:'../src', sourceFile:'kernel_nussionov.c', compiler:'llvm', nsteps:1000, language:'c', target:'host-intel', algo:'sa', metric:'performance', compilerFolderName:'/opt/clang+llvm-3.7.1-x86_64-linux-gnu-ubuntu-15.10/bin', seqlen:32, nexec:30, nr:-1, clean:1, passes:'', percent:2, append:''}" -nw -nci -b 2
```

### Using GCC

#### Minimum options (using defaults)
```
clava ../src/LaradLauncher.lara -av "('../src', 'kernel_nussionov.c', 'gcc')" -nw -nci -b 2
```
#### All options (with named parameters)
```
clava ../src/LaradLauncher.lara -av "{laradFoldername:'../src', sourceFile:'kernel_nussionov.c', compiler:'gcc', nsteps:1000, language:'c', target:'host-intel', algo:'sa', metric:'performance', compilerFolderName:'/opt/gcc-5.4.0/bin', seqlen:32, nexec:30, nr:-1, clean:1, passes:'', percent:2, append:''}" -nw -nci -b 2
```
