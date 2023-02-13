Introduction
============

**gstlearn** is a package developed by the Geostatistical Team of Mines
ParisTech. It contains several well established geostatistical functions
and methods and some prototypes based on ongoing research.

Loading RGeostats
=================

In each R session, you must load the **gstlearn** library

R used as a pocket calculator
=============================

    1+3

    ## [1] 4

    exp(log(sqrt(2^2)))

    ## [1] 2

Assignments
-----------

    a = 1
    b <- 3
    a + b

    ## [1] 4

    d = a + b
    d

    ## [1] 4

    d = 7
    d

    ## [1] 7

Working directory
-----------------

    getwd()

    ## [1] "/home/fors/Projets/gstlearn/gstlearn/doc/courses/r"

    ls()

    ## [1] "a" "b" "d"

    rm(d)
    ls()

    ## [1] "a" "b"

Quit
----

    q()

Save the workspace (y).

The content of your session is saved in the file **.Rdata** (with the
**.Rhistory**) of your working directory.

*Caution*: these files are hidden.

To continue, you must relaunch R (and load RGeostats again).

Classes in R
------------

    b <- 3
    class(b)

    ## [1] "numeric"

    col = "blue"
    class(col)

    ## [1] "character"

    reponse=FALSE
    class(reponse)

    ## [1] "logical"

Booleans
--------

    a = 1
    a < 2

    ## [1] TRUE

    comp = a < 2
    comp

    ## [1] TRUE

    comp & reponse

    ## [1] FALSE

    comp | reponse

    ## [1] TRUE

    !comp

    ## [1] FALSE

    a = 1
    a == 1

    ## [1] TRUE

    a != 1

    ## [1] FALSE

    a > 1

    ## [1] FALSE

    a >= 1

    ## [1] TRUE

Vectors
-------

    x = c(1,9,2,9,4,5)
    x

    ## [1] 1 9 2 9 4 5

    class(x)

    ## [1] "numeric"

    x[1]

    ## [1] 1

    x[c(1,3)]

    ## [1] 1 2

    x[-2]

    ## [1] 1 2 9 4 5

    2:5*3

    ## [1]  6  9 12 15

    x = c(1,9,2,9,4,5)
    x[2:5]

    ## [1] 9 2 9 4

    selection = x > 2
    selection

    ## [1] FALSE  TRUE FALSE  TRUE  TRUE  TRUE

    x[selection]

    ## [1] 9 9 4 5

Length of a vector

    length(x)

    ## [1] 6

Matrices
--------

    M = matrix(1:20,nrow=4,ncol=5)
    M[2,1]

    ## [1] 2

    M[2,1]=3
    M[2,1]

    ## [1] 3

    M[,1]

    ## [1] 1 3 3 4

    M[2,]

    ## [1]  3  6 10 14 18

    dim(M)

    ## [1] 4 5

Data frame
----------

    data = data.frame(M)
    class(data)

    ## [1] "data.frame"

    data[,3]

    ## [1]  9 10 11 12

    names(data)

    ## [1] "X1" "X2" "X3" "X4" "X5"

    names(data) = c("A","B","C","D","E")
    data

    ##   A B  C  D  E
    ## 1 1 5  9 13 17
    ## 2 3 6 10 14 18
    ## 3 3 7 11 15 19
    ## 4 4 8 12 16 20

    data

    ##   A B  C  D  E
    ## 1 1 5  9 13 17
    ## 2 3 6 10 14 18
    ## 3 3 7 11 15 19
    ## 4 4 8 12 16 20

    data$C

    ## [1]  9 10 11 12

    data[data$A>1,]

    ##   A B  C  D  E
    ## 2 3 6 10 14 18
    ## 3 3 7 11 15 19
    ## 4 4 8 12 16 20

Data Frame from a File
----------------------

Download the ASCII file called *Scotland\_Temperatures.csv* and store it
on your disk in the current working directory. For now on, the reference
files are stored on disk at an address referenced by the Environment
variable *GSTLEARN\_DATA*.

    getwd()
    filename = paste(Sys.getenv('GSTLEARN_DATA'),"Scotland",
                     "Scotland_Temperatures.csv",sep="/")
    temperatures = read.csv(filename,header=TRUE,na="MISS")
    class(temperatures)
    names(temperatures)

Some functions
--------------

    set.seed(123)
    rnorm(3)

    ## [1] -0.5604756 -0.2301775  1.5587083

    x=c(3,4,7,4,2)
    x

    ## [1] 3 4 7 4 2

    sort(x)

    ## [1] 2 3 4 4 7

    unique(x)

    ## [1] 3 4 7 2

    x=c(3,4,7,4,2)
    sum(x)

    ## [1] 20

    mean(x)

    ## [1] 4

    var(x)

    ## [1] 3.5

    sd(x)

    ## [1] 1.870829

Help of Functions
-----------------

Search for help on the function **image**.

    ?image

To get some help on a class :

    class?matrix

Some functions of interest
--------------------------

Function **outer**

    A = c(2,4,10)
    B = c(3,7)
    outer(A,B,"-")

    ##      [,1] [,2]
    ## [1,]   -1   -5
    ## [2,]    1   -3
    ## [3,]    7    3

    outer(A,B,"*")

    ##      [,1] [,2]
    ## [1,]    6   14
    ## [2,]   12   28
    ## [3,]   30   70

Loops
-----

Dimension and initialize a vector: **numeric**

Loop from 1 to 10: store loop index in a vector; then print the vector

    tab = numeric(10)
    for (i in 1:10)
    {
      tab[i] = i
    }
    print(tab)

    ##  [1]  1  2  3  4  5  6  7  8  9 10

Create your own function
------------------------

Generate a linear transform of the input argument

    my.func = function(x,a=2,b=4)
    {
      y = a * x + b
      y
    }

Play my function

    x = seq(0,20)
    y = my.func(x)
    plot(x,y,main="Plotting my function")

