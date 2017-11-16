# FlowCutter Code

This is the experimental code written during the development of FlowCutter. The code is provided as-is. It was never intended to be released and is thus poorly documented and will likely contain bugs. However, it might be of use to other researchers. If you use the code please cite:

* Graph Bisection with Pareto-Optimization
  Michael Hamann, Ben Strasser
  2016 Proceedings of the Eighteenth Workshop on Algorithm Engineering and Experiments (ALENEX) 

To get started run:

```bash
./build.py
./console interactive
```

Type `help` to get an overview of the available commands. Typing `details command` sometimes gives a detailed description of a command. If there is no description, you must read the source code.

The commands above should work on all Unix systems. On Windows, you will at least run into problems with directory separators.

License: The code in this repository is under BSD license. However, one can optionally link libraries, whose code is not copied in this repository, that have a GPL license. If you link these libraries, the code in this repository is also under GPL for the usage case. The relevant libraries are

* [[https://tiswww.case.edu/php/chet/readline/rltop.html|GNU readline library]] to enable tab-completion. By default it is is compiled in. You can define the macro NO_GPL to get rid of it. In this case, tab-completion does not work. Fortunately, all major functionality is untouched by this macro.

* [[http://algo2.iti.kit.edu/kahip/|KaHIP]]. By default, KaHip is not linked. You can build a variant with KaHIP support using `./build.sh`. This alternative build mode enables commands that compute nested dissection orders using KaHIP.

