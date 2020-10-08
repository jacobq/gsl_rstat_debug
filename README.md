# Runnings stat checks

Compares [`gsl`'s](https://www.gnu.org/software/gsl/) [running statistics](https://www.gnu.org/software/gsl/doc/html/rstat.html) mean/variance implementation with homemade one.
There is a [known numerical challenge](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance) regarding computing the variance.
See [Welford's online algorithm](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm)

```
sudo apt install libgsl-dev libgsl-dbg
```