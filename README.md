# STN-GPe-synaptic-rewiring

## Rhythmic modulation of subthalamo-pallidal interactions depends on synaptic rewiring through inhibitory plasticity

**Abstract**: ‎Rhythmic stimulation offers a paradigm to modulate brain oscillations and‎, ‎therefore‎, ‎influence brain function‎. ‎A huge body of evidence indicates that reciprocal interactions between the neurons of the subthalamic nucleus (STN) and globus pallidus externus (GPe) play a central role in the emergence of abnormal synchronous beta (15-30 Hz) oscillations in Parkinson's disease (PD)‎. ‎Proliferation of the inhibitory GPe-to-STN synapses following dopamine loss adds to this pathological activity‎. ‎Rhythmic modulation of the STN and/or GPe‎, ‎e.g.‎, ‎by deep brain stimulation (DBS)‎, ‎can restore physiological patterns of activity and connectivity‎. ‎Here‎, ‎we tested whether dual STN-GPe targeting by rhythmic stimulation can modulate pathologically strong GPe-to-STN synapses through inhibitory spike-timing-dependent plasticity (iSTDP)‎. ‎More specifically‎, ‎we examined how time-shifted paired stimuli delivered to the STN and GPe can lead to inter-population synaptic rewiring‎. ‎To that end‎, ‎we first theoretically analysed the optimal range of stimulation time shift and frequency for effective synaptic rewiring‎. ‎Then‎, ‎as a minimal model for generating subthalamo-pallidal oscillations in normal and PD conditions we considered a biologically inspired STN-GPe loop comprised of conductance-based spiking neurons‎. ‎Consistent with the theoretical predictions‎, ‎rhythmic stimulation with appropriate time shift and frequency modified GPe-to-STN interactions through iSTDP‎, ‎i.e.‎, ‎by long-lasting rewiring of pathologically strong synaptic connectivity‎. ‎This ultimately caused desynchronising after-effects within each population such that excessively synchronous beta activity in the PD state was suppressed‎, ‎resulting in a decoupling of the STN-GPe network and restoration of normal dynamics in the model‎. Decoupling effects of the dual STN-GPe stimulation can be realized by time-shifted continuous and intermittent stimuli as well as monopolar and bipolar simulation waveforms. ‎‎‎Our findings underlines the role of plasticity in shaping long-lasting stimulation effects and may contribute to the optimisation of a variety of multisite stimulation paradigms aimed at a reshaping of diseased brain networks by targeting plasticity.

## Usage

- These codes reproduce data used to generate the figures in the manuscript. ```STN_GPe_network_main.cpp``` reproduces the key results, ```iSTDP_profile.cpp``` reproduces the symmetric iSTDP learning window, ```stimulation_current.cpp``` reproduces monopolar-continuous, biopolar-continuous and monopolar-burst stimulation waveforms, ```structural_connectivity.cpp``` reproduces the structural connectivity of the STN-GPe network with random connections, and ```theoretical_predictions.cpp``` reproduces the theoretical predictions of the iSTDP-induced synaptic changes in the model.

#### Compilation

- Compile and execute either of the scripts with the following syntax:

```
g++ -std=c++11 file_name.cpp
./a.out
```

#### Requirements

- C++11,
- g++ compiler,
- ```random_number.cpp``` and ```spike_statistics.cpp```.

## Citation

- Madadi Asl, M., & ‎Lea-Carnall, C. A. (2025). Rhythmic modulation of subthalamo-pallidal interactions depends on synaptic rewiring through inhibitory plasticity. Physical Review Research, 7: 023128. [https://doi.org/10.1371/journal.pcbi.1010853](https://doi.org/10.1371/journal.pcbi.1010853).
