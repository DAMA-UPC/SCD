SCD
===

This program is an implementation of the community detection algorithm described in the paper titled

[High quality, scalable and parallel community detection for large real graphs.](http://www.dama.upc.edu/publications) Arnau Prat-Pérez, David Dominguez-Sal, Josep-Lluis Larriba-Pey - WWW 2014.


Compile
===

SCD uses CMake 2.8.2 or greater to compile. In order to build SCD, move to SCD directory and type:

```
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
``` 

This will create a build directory into the SCD folder tree, and configure and build SCD in Release mode.
In order to compile SCD in Debug mode, please replace the last two lines of the snippet above by:

```
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
``` 

Execution
===

To execute SCD type, move to the build folder and type:

```
./scd -f [network file name]
```

where the [network file name] contains the network with an edge per line, and each edge is represented as a pair of numeric identifiers. 
IMPORTANT: each edge is interpreted as an undirected edge and can only appear once. 
For example, if "1 2" appears in the file, then "2 1" cannot appear too. The next snipped shows a valid network file:

```
1 2
4 5
3 4
1 3
```

As an example, type:

```
./scd -f ./network.dat
```

This will run the program, and will output the communities found at network.dat to "./communities.dat", which contains
a community per line, represented as a list of identifiers. network.dat contains a network composed by two cliques of size 5 linked by a single edge. Therefore, communities.dat outputs:

``` 
1 2 3 4 5
6 7 8 9 10
```

Now, we summarize the options of the program:

  *  -f [netork file name] : Specifies the network file.
  *  -o [output file name] : Specifies the output file name.
  *  -n [number of threads]: Specifies the number of threads to run the algorithm.





