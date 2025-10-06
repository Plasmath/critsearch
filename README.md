# critsearch

Critsearch is a shiny new command-line tool designed to find new polytopes (such as [uniform polytopes](https://polytope.miraheze.org/wiki/Uniform_polytope) and [noble polytopes](https://polytope.miraheze.org/wiki/Noble_polytope)) by searching for [critical points](https://polytope.miraheze.org/wiki/Critical_point) within a given orbit type. This project is written in C++ and contains the [nlohmann/json](https://github.com/nlohmann/json) repository as a dependency in order to read JSON files.

An 'orbit type', informally, is a collection of orbits of a symmetry group that can be continuously morphed between each other.

Currently, the program is in its initial stages of development, so it contains many limitations and bugs (apologies for my terrible code!). Contributions to the library and/or codebase are welcome.

## Running the program

The program can be compiled with g++ by using the following command. This command should work on any Unix-like environment - I ran this on Windows through the [Cygwin](https://cygwin.com/) terminal.

`g++ critsearch.cpp -I ./ -std=c++11 -o critsearch`

One may then simply run the compiled `critsearch` executable. Once in the program, use the `help` command to see a list of commands.

## Current algorithms

Critsearch currently supports the following algorithms:

- `qunivar` - short for 'quadratic univariate search'. As the name suggests, it searches for uniform/scaliform polytopes while only using the quadratic formula to find critical points, and so it does not work when there is more than one free variable.

## Format for orbit types

Critsearch stores its search spaces (orbit types) with JSON files that describe

a. The symmetry group of the orbits
a. The initial points of the orbits that, along with the symmetry group, generate the full set of points within an orbit
a. OFF file data that allows critsearch to automatically export its results into the OFF format.

This is stored through a JSON object with the following keys:

- `name` : The name of the orbit type as a string (e.g. truncated cube).
- `author` (optional) : Name of the JSON's author.
- `dim` : The dimension of the orbit as an integer.
- `num-variables` : An integer indicating the number of free variables (a.k.a. degrees of freedom) in the orbit type (NOTE: No algorithms support more than one free variable at the moment).
- `node-activations` : An array of strings, each indicating a variable symbolic expression (e.g. `"0"` or `"-sqrt(2)/2"` or `"x"`) of how 'active' each generating vector should be (see `generating-vectors`).
- `generating-matrices` : Generators for the matrix representation of the symmetry group. Each entry in a matrix is a string that could be a constant symbolic expression (e.g. `"0"` or `"-sqrt(2)/2"`, but not `"x"`).
- `generating-vectors` : Tells critsearch which vectors count as 'activating a node' (each entry is also a constant symbolic expression). When the program stars it will multiply each vector with its corresponding variable in `node-activations`.
- `off-header` : When exporting, put this string before the vertex data (e.g. `4OFF\n288 336 576 48`). Note the usage of the newline escape sequence `\n`.
- `off-body` : When exporting, put this string after the vertex data (a.k.a. the data regarding faces, cells, tera, etc.). This usually contains the bulk of the data in the JSON file.

Examples can be found in the library folder, which is where the program will look for JSON files when importing. Currently a 'constant expression' consists only of fractions & square roots.