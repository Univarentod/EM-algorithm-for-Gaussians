# EM-algorithm-for-Gaussians
Implementation of EM algorithm for mixture of Gaussians

## How to Use

These notes would focus on the input and output of the EM function. The "black box".

### Input

There are 2 parameters.

1. input_data

The format of input data is limited as a list of lists. This means if you want to use a set of data as input, say, 1, 2, 3, 4, 5, you'll have to pack them up into a list before feed them to the function:

```
input_data = [[1], [2], [3], [4], [5]]
```

2. k

This parameter determines how many clusters will the function cluster your input data into. Assign a positive integar to this parameter. Or you can just ignore it, it would be assigned to 2 by default.

### Output

The function would return the clustered input data and the most possible Gaussian distributions for these data. You can use `output[0]` to get the clustered input data, and `output[1]` to get the information of the clusters. However, when k == 1 you will only get the input data and when k < 1 you'll get nothing.

## Accuracy

When your input data are symmetrically distributed, the clusters may all converge to some centroids. When your input data are too concentrated, its accuracy would drops rapidly. All of these situations would make the function perform no better than randomly clustering the input data.

## Potential Problems when Reading my Codes

You may find some codes illogical. Hope these would be able to answer your question.

1. Initialization

To calculate the pdf of multivariate normal distribution, we would need a positive definite covariance matrix. Therefore, it would be crucial for us to build up such a matrix during the initialization process. My solution is:

```
i.   randomly pick points from the input data to determine expectations of the clusters;
ii.  randomly pick a point from the input data to form a covariance matrix, using the formula of covariance matrix(check the one on Wikipedia yourself);
iii. collect positive eigenvalues from the matrix
iv.  repeat ii to iii until there are enough eigenvalues to form a corresponding positive definite matrix
```

2. E-step

Usually the pdf would come up with a very very small number that the computer would make it 0 due to the precision loss. So in the computation of weights of ownership, I would do the computation in log space. Also, I would divide each element with the first element so that they would get values relatively close to 0(in log space), which makes it possible to add the numbers(in real space).

## Acknowledgments

* [The Standford's CS229 classnotes 8](http://cs229.stanford.edu/notes/cs229-notes8.pdf) explains how the EM works
* [Second partial derivative test](https://www.khanacademy.org/math/multivariable-calculus/applications-of-multivariable-derivatives/optimizing-multivariable-functions/a/second-partial-derivative-test) is useful to help understand the classnote above

## Contact Me

If you have any problems, welcome to send e-mails to [my mailbox](mailto:657711905@qq.com)
