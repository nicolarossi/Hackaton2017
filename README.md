
Note:
Sorry to my english mistakes, to correction/question/improvement and whatever you want feel free to contact me on gmail.com using name and surname .

# Introduction

In this document are showed the test results about "Black-Sholes" challenge on The DutchCppGrp - Hackaton 2017 *Thanks to organizator for the beautiful experience.

The result are provided using avg and sigma values calculated on 100 samples present in the file "Sample_SPX_20151001_to_20151030_WITH_SPXW.csv"
which  only the "^SPXW," grepping rows present in the original file provided us.

# Exercise 1

Optimize the call_price function and the uploading 

To test  different version of the code, you can play to define in the source code or compile-time this macro

*  FUNCTION
*  USING_MACRO
*  REORDER_1
*  REORDER_2
*  INDIRECT_INDEX
*  USING_REDUCTION

## FUNCTION

This was the original version of the code .

By example:
```
double d1(double S, double K, double r, double v, double T)
{
  return (log(S / K) + r*T) / (v*(pow(T, 0.5))) + 0.5*v*(pow(T, 0.5));
}

```

## USING_MACRO

This was create substituting the function-like definition with a macro-like definition to "try" to
help the compiler to optimize that calculum

```

#define norm_cdf( value) ( 0.5 * erfc(-value * M_SQRT1_2))
#define d1( S,  K,  r,  v,  T) ( (log(S / K) + r*T) / (v*(pow(T, 0.5))) + 0.5*v*(pow(T, 0.5)))
#define d2( S,  K,  r,  v, T)  ( d1(S, K, r, v, T) - v*(pow(T, 0.5)))

```

## REORDER_1

Looking the assembly generated we see that are present multiple call at pow function so we try to
force the compiler to create better code.

```

#define d1( S,  K,  r,  v,  T,v_pow) (( (log(S / K) + r*T) / (v_pow)) + 0.5*v_pow )

```

## REORDER_2

Same idea of last step but observing that S,K,r, T are constant in our research.



```

#define norm_cdf( value) ( 0.5 * erfc(-value * M_SQRT1_2))
#define d1( S,  K,  r,  v,  T,v_pow,log_S_K_p_rT) (( (log_S_K_p_rT) / (v_pow)) + 0.5*v_pow )
#define d2( S,  K,  r,  v, T,v_pow,log_S_K_p_rT)  (( (log_S_K_p_rT) / (v_pow)) - 0.5*v_pow )

```

 the resulting code will be:

```

//--- Outside loop on stock rows ...
exp_rT=exp(r*T);
c_pow=pow(T,0.5);

for_each_row_in_input_file

//--- ... inside the loop on stock row
double v_pow=v*c_pow;

log_S_K_p_RT=log(S/K)+r*T;
best_finded=find_volatility(call_real, S, K,  r,  T,  a,b, dv);
    
...

double _call_price(double S, double K, double r, double v, double T) 
  return S * norm_cdf(d1(S, K, r, v, T,v_pow,log_S_K_p_RT)) - K*exp_rT * norm_cdf(d2(S, K, r, v, T,v_pow,log_S_K_p_RT));


```

# Exercise 2


Explore the domain of volatility to find the value of "v" that generate a "call_price" .


# dv static

In this version of the software we search volatility using a costant step loop to find a volatility "v" so the value call_price(S,K,v,T) is satisfied.

# dv dynamic

Using a $(sorry_i_have_forgotten_the_name)'s idea, we use a dynamic size of "dv"  but limited to the minimum value.

The ratio v/dv could be constant, this implies that little value of v  implies dv = 1e-4 ; bigger value of v implies dv_new = v*dv
so the values of dv will be (by example)

| Volatility | delta to search |
|---|---|
| 0 | 0.001 |
| 0.001 | 0.001 |
| 0.002 | 0.001 |
| ... | ... |
| 10.01 | 0.01 |
| 10.02 | 0.01 |
| 10.03 | 0.01 |
| ... | ... |
| 100.1 | 0.1 |
| 100.2 | 0.1 |
| 100.3 | 0.1 |



# Indirect index

To enable a correct and easy parallelization we
1) create the vector of "dv" step and the starting point of the search
2) parallelize the work with a simple division of that vector.
3) gather the threads local solutions.

# No reduction

To gather the solution of the 2-th step of the parallelization stage,
we can create a vector on which every thread wrote its local solution
of the problem and after the parallel section the Return Value vector
is explored.


```

 int who_i_am=omp_get_thread_num();
    //    std::cerr<<" debug "<<i<<"\t"<< tmp_rv.call<<"\t\t";
    if (tmp_rv<rv[who_i_am]) {
      rv[who_i_am]=tmp_rv;
      //      std::cerr<<" cambiamo ";
    }
  ....
    
  for (int i=0;i<omp_get_max_threads();i++){
    if (rv[i]<rv[index_min]){
      index_min=i;
    }
  }

  return rv[index_min];
  
```

# Using reduction 

For every thread the element in the rv vector of its ownership is adjacent to the rv element of
the neighbor thread.

* rv[0] is the rv-element for the thread 1
* rv[1] is the rv-element for the thread 2
* rv[2] is the rv-element for the thread 3

Due to the spacial locality, in the memory hierarchy system, more element of an array can be stored
in the same row of cache.

This is an advantage in the reading stage but in the writing stage can become a disadvantage, due to
the necessity to synchronize the data between cache and memory propagating the write operation to
the lower level of the hierarchy memory.

So it can be more powerful to reduce the number of writes in a vector shared between more threads,
especially when the writing is made on the same cache block.

To avoid this sharing we used the OpenMP reduction clause 

```

 solution_t rv_min;

 #pragma omp declare reduction ( best_solution : solution_t : omp_out = omp_in < omp_out ? omp_in : omp_out ) 
 #pragma omp parallel for private(tmp_rv,v,dv) reduction(best_solution : rv_min)

 for (int i=0;i<step.size();i++) {
    v=(step[i].v);
    tmp_rv.call_price(S, K,  r,  v,  T);

    if (tmp_rv<rv_min) {
      rv_min=tmp_rv;
    }
    dv=step[i].dv;
  }
  
```

The test result in this case are very surprising and unexpected.


# Test result

Due to the nature of the domain it would have been more efficient, easier but very less fun to parallelize the problem assigning every rows in the file to a different thread.

To personal challenge I tried to parallelize every single stock and obviously are less efficient, but more interesting.


## on i5 5500U 

| Optimization   | avg (msec/stock) | sigma  | Note   |
|-------|---|---|---|---|
|  dv Static  | 229.34  | 2.22  |   |   |
|  dv Dynamic  | 7.67  | 0.77  |  INITIAL_DV=1e-4 |   |
|  Function  | 7.64  | 1.08  |   |   |
|  Using_Macro | 7.59  | 0.76  |   |   |
|  Reorder_1 |  7.33 | 0.49  |   |   |
|  Reorder_2 |  3.41 |  0.64 |   |   |


## Scalability

### on i5 5500U 

Using dv_dynamic + Reorder_2

For every cell we report avg &plusmn; sigma valuem and efficiency ( time with 1 threads / ( n * time_n threads) )


|  Test | 1 thread | | 2 threads | | 4 threads |  |8 threads | |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|  With INDIRECT_INDEX | 3.29 &plusmn; 0.88  |  n/a   | 2.05 &plusmn; 0.21  | 80%  | 2.35 &plusmn; 1.00  | 35%  | 1.56 &plusmn; 0.84  | 26% |
|  with reduction clause | 3.38 &plusmn; 0.86  | n/a  | 2.16 &plusmn; 0.65  | 78%  | 2.51 &plusmn; 1.19 |  33% |  1.54 &plusmn; 0.81 | 27%  |


### on i5  XeonPhi

|  Test | 1 thread | | 2 threads | | 4 threads |  |8 threads | |16 threads | |32 threads | |64 threads | |128 threads | |
|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
|  With INDIRECT_INDEX | 79.84 &plusmn; 2.29  |  n/a   | 42.12 &plusmn; 1.29  | 93%  | 23.58 &plusmn; 0.66  | 84%  | 14.17 &plusmn; 0.9  | 69% | 7.11 &plusmn; 1.09  | 69.44% | 3.19 &plusmn; 0.41 ** | 77% | 1.26 &plusmn; 0.81 **  | 97% |1.08 &plusmn; 0.63 ***  | 57% |
|  with reduction clause | 80.3 &plusmn; 2.60  | n/a | 42.06 &plusmn; 0.70  | 95% | 23.13 &plusmn; 1.01  | 87% |14.14 &plusmn; 0.64  | 71% |7.1 &plusmn; 0.9  | 71% |3.019 &plusmn; 0.4 **  | 83% |1.26 &plusmn; 0.75 ** | 98%  |1.016 &plusmn; 0.61 ** | 62%  |


** the measure is made with 1000 samples
*** the measure is made with 10000 samples





