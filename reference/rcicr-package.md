# rcicr: Reverse-Correlation Image-Classification Toolbox

Generate stimuli and analyze data of reverse correlation image
classification experiments (psychophysical tasks aimed at visualizing
cognitive mental representations of faces). For the method see Dotsch
and Todorov (2012)
[doi:10.1177/1948550611430272](https://doi.org/10.1177/1948550611430272)
; for a practical primer see Brinkman, Todorov and Dotsch (2017)
[doi:10.1080/10463283.2017.1381469](https://doi.org/10.1080/10463283.2017.1381469)
.

## Details

[`vignette("reverse-correlation-walkthrough")`](https://rdotsch.github.io/rcicr/articles/reverse-correlation-walkthrough.md)
works through a complete experiment end to end. In brief:
[`generateStimuli2IFC`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
creates the stimuli together with the `.Rdata` file recording how they
were made, and
[`generateCI2IFC`](https://rdotsch.github.io/rcicr/reference/generateCI2IFC.md)
or
[`batchGenerateCI`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
turns participants' responses into classification images.
`citation("rcicr")` prints the reference to cite.

## References

Dotsch, R., & Todorov, A. (2012). Reverse correlating social face
perception. *Social Psychological and Personality Science*, *3*(5),
562-571.
[doi:10.1177/1948550611430272](https://doi.org/10.1177/1948550611430272)

Brinkman, L., Todorov, A., & Dotsch, R. (2017). Visualising mental
representations: A primer on noise-based reverse correlation in social
psychology. *European Review of Social Psychology*, *28*(1), 333-361.
[doi:10.1080/10463283.2017.1381469](https://doi.org/10.1080/10463283.2017.1381469)

Dotsch, R., Wigboldus, D. H. J., Langner, O., & Van Knippenberg, A.
(2008). Ethnic out-group faces are biased in the prejudiced mind.
*Psychological Science*, *19*(10), 978-980.
[doi:10.1111/j.1467-9280.2008.02186.x](https://doi.org/10.1111/j.1467-9280.2008.02186.x)

## See also

Useful links:

- <https://rdotsch.github.io/rcicr/>

- <https://github.com/rdotsch/rcicr>

- Report bugs at <https://github.com/rdotsch/rcicr/issues>

## Author

**Maintainer**: Ron Dotsch <rdotsch@gmail.com>
