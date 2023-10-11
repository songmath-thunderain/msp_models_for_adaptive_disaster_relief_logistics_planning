# [Multistage Stochastic Programming with a Random Number of Stages: Applications in Hurricane Disaster Relief Logistics Planning](https://arxiv.org/pdf/2201.10678.pdf)
## authors:
  1. Murwan Siddig: RWTH Aachen University, [siddig@dpo.rwth-aachen.de](siddig@dpo.rwth-aachen.de)
  2. Yongjia Song: Clemson University, [yongjis@clemson.edu](yongjis@clemson.edu)
## abstract:
We consider a logistics planning problem of prepositioning relief items in preparation for an impending hurricane landfall. We model the problem as a multi-period network flow problem where the objective is to minimize the total expected logistics cost of operating the network to meet the demand for relief items. We assume that the hurricane’s attributes evolve over time according to a Markov chain model, and the demand for relief items is calculated based on the hurricane’s attributes (intensity and location) at the time of landfall, if one is to be made at all. We introduce a fully adaptive multi-stage stochastic programming (MSP) model that allows the decision-maker to adapt their logistics decisions over time according to the evolution of the hurricane’s attributes. In particular, we develop a novel extension of the standard MSP model to address the challenge of having a random number of stages in the planning horizon due to the uncertain landfall or dissipation time of the hurricane. Our work is an extension of the approach proposed in Guigues
(2021) to the setting where the underlying stochastic process is modeled as a Markov chain. We benchmark the performance of the adaptive decision policy given by the MSP models with alternative decision policies, such as the ones obtained from the static and rolling-horizon two-stage stochastic programming approaches. Our numerical results provide key insights into the value of MSP in the hurricane disaster relief logistics planning problem.

## paper: preprint version [link](https://arxiv.org/pdf/2201.10678.pdf)

## Example run command: 
python main.py -p solveParams.yaml -d 1 -a 1 -k 57 -o 1000 -ni 3 -nj 10 -t 0.05 -s 1

## Get help on all different arguments to use:
python main.py -h