import scipy
import scipy.io
import sys
import scikits.sparse.cholmod

filenames=[sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10]]
A=scipy.io.mmread(filenames[0])
n=A.shape[0]
offtree=scipy.sparse.triu(A,2)
U=scipy.sparse.triu(A,1)
tree=U-offtree
tree=tree+tree.transpose()
diag=tree.sum(1)
diag=-diag
tree=tree.tolil()
tree.setdiag(diag)

scipy.io.mmwrite(filenames[1],tree,precision=9)

nnz=U.nnz
perm=scipy.argsort(-offtree.data)
index=scipy.arange(n-1)
fraction=[.01,.1,.15,.2]
for i in xrange(0,4):
    extra=scipy.floor(fraction[i]*nnz)
    tree=U-offtree
    subdata=offtree.data[perm][0:extra]
    subrow=offtree.row[perm][0:extra]
    subcol=offtree.col[perm][0:extra]
    treeplus=scipy.sparse.coo_matrix((subdata,(subrow, subcol)), (n,n))
    treeplus=treeplus+tree
    treeplus=treeplus+treeplus.transpose()
    diag=-treeplus.sum(1)
    treeplus=treeplus.tolil()
    treeplus.setdiag(diag)
    nonsingular=treeplus.tocsr()[index,:].tocsc()[:,index]
    scipy.io.mmwrite(filenames[i+2],nonsingular,precision=9)
    factor=scikits.sparse.cholmod.cholesky(nonsingular)
    L=factor.L()
    scipy.io.mmwrite(filenames[i+6],L,precision=9)
    
