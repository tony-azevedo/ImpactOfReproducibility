Theory project - Functional impact of the reproducibility of sparse signals in a convergent, thresholding network
Purpose:  To describe the fidelity of quantal responses, and to examine how this fidelity limits the output of a network that pools across thresholding elements before thresholding the output.

Hypothesis: 
Sequential non-linearities that act on pooled, sparse signals are positioned in such a way as to balance preservation of rare but relevant signals at the cost of rejecting signals (actual events) that may be confused with noise.  Altering the signal distribution in a way that changes the overlap between signals and noise will disrupt the performance of this system.

Methods:

   * Use single photon responses measured in mouse rods.
   * Use flash responses in ganglion cells at visual threshold to constrain the output of the network.
   * Build a network that features two thresholding non-linearites, one that acts on the sparse input, and another that acts on the pooled output of that non-linearity.
   * Simulate inputs drawn from the real distributions and from modified distributions

Model parameters:

   * Signal and noise distributions in rod responses (sigma_D, sigma_i, A_bar)
   * Thresholding
   * Linear transfer functions and pooling
   * Thresholding nonlinearity

Analysis:

   * Measure the performance of the system in a dim flash detection task
   * Measure the discriminability of flashes of different intensities
Things to examine: 

   * How would a linear system respond to changes in input signal variability?  Our initial idea is that the integral of the response (charge) and it's associated fidelity would determine the system response.
   * How do spectral differences in signal and noise affect their distributions?  Do we match filter the responses right away?  What if we use the entire response? How does this affect my data?
   * Can reproducibility can be represented as the spread in the signal amplitude if the appropriate match filter is applied?
   * Can a nonlinearity remove variability present in the initial signal?
Caveats:


   * It will be tough to model all of the elements in this circuit.  

