# ❀ Swarm Intelligence in Bioinformatics ❀
#### •● Table of Content
- Traditionals
  - [Particle Swarm Optimization](#-particle-swarm-optimization-pso)
  - [Ant Colony Optimization](#-ant-colony-optimization-aco)
  - [Artificial Bee Colony](#-artificial-bee-colony-abc)
- New Additions
  - [Grey Wolf Optimization](#-grey-wolf-optimization-gwo)
  - [Whale Optimization Algorithm](#-whale-optimization-algorithm-woa)
  - [Social Spider Optimization](#-social-spider-optimization-sso)
  - [Spider Monkey Optimization](#-spider-monkey-optimization-smo)
  - [Moth Search Algorithm](#-moth-search-algorithm-msa)
  - [Monarch Butterfly Optimization](#-monarch-butterfly-optimization-mbo)

<br><br>

#### •● Particle Swarm Optimization (PSO)
```
Algorithm 1: Particle Swarm Optimization

```
<br><br>

#### •● Ant Colony Optimization (ACO)
```
Algorithm 2: Ant Colony Optimization

```
<br><br>

#### •● Artificial Bee Colony (ABC)
Artificial Bee Colony was first proposed by Karaboga in 2005. It is an optimization algorithm inspired by the intelligent foraging behavior of a honey bee swarm for finding an optimal solution. Three essential components of forage selection are:
1)	food sources: places with a high amount of nectar 
2)	employed foragers: employed bees exploit the food source, calculate the nector amount, and carry the optimal path information back to hive.
3)	unemployed foragers: unemployed bees consist of two groups of bees – onlookers and scouts. Onlookers wait in the hive for the information that is shared by the employed bees. Scout bees are translated from employed bees and always continues seeking new food sources near the hive. <br>
Comparing Artificial Bee Colony algorithm with Particle Swarm Optimization and Ant Colony Optimization methods, it carries the advantages of fast convergence, high flexibility, and strong robustness. 
```
Algorithm 3: Artificial Bee Colony
1:	Initialize food sources
2:	while termination criteria is not met do
3:	  for each employed bee
4:	    Produce new solution
5:	    Calculate the fitness value
6:	    Apply greedy selection 
7:	    Calculate the probability value
8:	  end
9:	  for each onlooker bee
10:	    Select a solution 
11:	    Produce new solution
12:	    Calculate the fitness value
13:	    Apply greedy selection 
14:	  end
15:	  if an abandoned solution for the scout exists then
16:	    Replace it with a new solution at random
17:	  end
18:	  Register the best solution
19:	end
```
<br><br>

#### •● Grey Wolf Optimization (GWO)
Reference: Grey Wolf Optimizer, Advances in Engineering Software, Volume 69, March 2014, Pages 46-61, DOI: http://dx.doi.org/10.1016/j.advengsoft.2013.12.007 
```
Algorithm 4: Grey Wolf Optimization Algorithm
Initialize the grey wolves population Xi(i = 1; 2; :::; n)
Initialize a, A and C
1: Calculate the fitness of each search agent
2: Xa = the best search agent
3: Xb = the second best search agent
4: Xc = the third best search agent
5: while t < Max number of iterations do
6:    for each search agent do
7:      Update the position of the current search agent by Eq. (1)
8:    end 
9:    Update a, A, and C
10:   Calculate the fitness of all search agents
11:   Update Xa, Xa and Xc
12:   t = t + 1
13: end 
14: return Xa
```
<br><br>

#### •● Whale Optimization Algorithm (WOA)
Reference: The Whale Optimization Algorithm, Advances in Engineering Software, in press, 2016. DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008 
```
Algorithm 5: Whale Optimization Algorithm

```
<br><br>

#### •● Social Spider Optimization (SSO)
Reference: E. Cuevas et al. Swarm optimization algorithm inspired in the behavior of the social-spider, Expert Systems with Applications, 40 (16), (2013), pp. 6374-6384. DOI: http://arxiv.org/abs/1406.3282
```
Algorithm 6: Social Spider Optimization

```
<br><br>

#### •● Spider Monkey Optimization (SMO) 
Spider Monkey Optimization is inspired by fission-fusion social foraging behavior of Spider Monkeys. It can be applied in various areas such as remote sensing image classification, wireless protocol design, and electromagnetics and antenna optimization. The SMO algorithm follows the properties of self-organization and division of labor for obtainning intelligent swarming behaviors. Spider monkeys determine their positions by learning from self experience, local leaders and a global leader. It provides positive feedback of self-organization. Members of each group are spread to multiple directions for seeking a food source. If the global leader gets trapped, the group will divide into smaller subgroups and continue repeating the previous steps. Local and global leaders’ update their decisions according to negative feedback from the local leader limit and global leader limit. This process is known as a division of labor property. The SMO algorithm consists of 6 phases:
1)	Local Leader phase (LLP)
2)	Global Leader phase (GLP)
3)	Local Leader Learning phase (LLL)
4)	Global Leader Learning phase (GLL)
5)	Local Leader Decision phase (LLD)
6)	Global Leader Decision phase (GLD)<br>

```
Algorithm 7: Spider Monkey Optimization

```
<br><br>

#### •● Moth Search Algorithm (MSA)
Reference: G. Wang. "Moth search algorithm: a bio-inspired metaheuristic algorithm for global optimization problems." Memetic Computing (2016): 1-14. 
```
Algorithm 8: Moth Search Algorithm 

```
<br><br>


#### •● Monarch Butterfly Optimization (MBO)
Reference: G. Wang, X. Zhao, and S. Deb. "A novel monarch butterfly optimization with greedy strategy and self-adaptive." In Soft Computing and Machine Intelligence (ISCMI), 2015 Second International Conference on, pp. 45-50. IEEE, 2015.
```
Algorithm 9: Monarch Butterfly Optimization

```
<br><br>


