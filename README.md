<h1 align="center">
	<strong>
	MaxRank
	</strong>
	<br>
</h1>

Original MaxRank implementation in C++, converted to CMake and patched up. I would like to state that **THIS ISN'T MY CODE**.

## Requisites
The program calls two external binaries that aren't built by CMake since they come from the [qHull](https://github.com/qhull/qhull) repository. These two binaries are:
* qhalf
* qconvex

located in the *bin* folder. I've left my own binaries in as an indication, if they don't work on your system you must build and replace them.

## Usage
Run the MaxRank computation by calling:
```console
myQuadTree.exe -p 4096 -d DIM -i tmp/idx -f path/to/datafile.txt -q path/to/queryfile.txt -r N -m BA -h 7 -t 0 -o 0 -v
```
* **DIM** must be set to the dimensionality of the data.
* **datafile.csv** is the file containing all data points. It must be formatted in a weird way, check the iPython notebook.
* **queryfile.csv** is the file contaning the "IDs" of the points to compute the MaxRank of.
* **N** must be set to the number of points listed in the queryfile.

What the other option means can be seen by running the executable without parameters.

Example run:
```console
myQuadTree -p 4096 -d 3 -i tmp/idx -f examples/Test3D50/data.txt -q examples/Test3D50/query.txt -r 50 -m BA -h 7 -t 0 -o 0 -v
```

## Output
The output of the computation is written in the **myout.txt** file, which is located in the root project folder, containing the computed MaxRank of each point listed in the queryfile. The MaxRank is counted starting from 0 instead of 1.

## Notes
* The original, untouched implementation can be found in the **maxrank.rar.old** archive.
* Do not delete the **tmp** folder: the program won't run if you do.
* Regarding the *-m* paramenter, only the "BA" option actually works.
* Similarly for the *-o* parameter, the "optimized node intersection" doesn't give any meaningful results.
* The program considers the highest values as the best ones (the optimal point is (1, ..., 1)) which is the opposite w.r.t. the Python implementation. Consider flipping the values of your data points.
* The iPython notebook **Prep_MaxRank.ipynb** will format your data in the way the program requires it to be. The formatting paradigm is the following:
    * Indices must start from 1 and must be continuous.
    * Points are actually seen as volumes, presumably due to the usual implementation of R-trees, so they must be formatted as: id, x, y, z, x, y, z (along with some padding).
    * The maximum precision allowed is 10^-6

## Author
*Original author is unknown*