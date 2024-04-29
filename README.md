Accelerated survival analysis
================

Os modelos de regressão de sobrevivência de classe acelerada são um
grupo de técnicas estatísticas usadas para analisar dados em que a
variável de resultado é o tempo até um evento, também conhecidos como
dados de sobrevivência. Estes modelos diferem dos modelos de
sobrevivência padrão por incorporarem o efeito das covariáveis na
própria escala de tempo, em vez de apenas na taxa de risco. Aqui está
uma análise dos três tipos principais (AFT, AH, AO):

### 1 - Modelo de tempo de falha acelerado (AFT model)

O *AFT model* é uma abordagem paramétrica (ou semiparamétrica) que
oferece uma opção robusta de regressão em análise de sobrevivência. A
ideia central por trás dos modelos AFT é que as covariáveis atuam
multiplicativamente no log do tempo até a ocorrência de algum evento de
interesse, acelerando ou desacelerando o processo de falha (KALBFLEISCH;
PRENTICE, 2002). Assim, para a
![i](https://latex.codecogs.com/svg.image?i "i")-ésima observação, um
*AFT model* é caracterizado por sugerir que
<h1 align="center">

![\log(T_i) = \mathbf{x}\_i^{\top}\pmb{\beta} + \nu_i, \\\\\\\\i=1, \dots, n,](https://latex.codecogs.com/svg.image?%5Clog%28T_i%29%20%3D%20%5Cmathbf%7Bx%7D_i%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%20%2B%20%5Cnu_i%2C%20%5C%20%5C%20%5C%20%5C%20i%3D1%2C%20%5Cdots%2C%20n%2C "\log(T_i) = \mathbf{x}_i^{\top}\pmb{\beta} + \nu_i, \ \ \ \ i=1, \dots, n,")

</h1>

em que

- ![T_i](https://latex.codecogs.com/svg.image?T_i "T_i"):
  ![i](https://latex.codecogs.com/svg.image?i "i")-ésima observação da
  variável aleatória que define o tempo de falha
  (![T\>0](https://latex.codecogs.com/svg.image?T%3E0 "T>0"));
- ![\mathbf{x}=(x_1, \dots, x_k)^{\top}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bx%7D%3D%28x_1%2C%20%5Cdots%2C%20x_k%29%5E%7B%5Ctop%7D "\mathbf{x}=(x_1, \dots, x_k)^{\top}"):
  Vetor de covariáveis com dimensão
  ![k](https://latex.codecogs.com/svg.image?k "k") para a
  ![i](https://latex.codecogs.com/svg.image?i "i")-ésima observação;
- ![\pmb{\beta}=(\beta_1, \dots, \beta_k)^{\top}](https://latex.codecogs.com/svg.image?%5Cpmb%7B%5Cbeta%7D%3D%28%5Cbeta_1%2C%20%5Cdots%2C%20%5Cbeta_k%29%5E%7B%5Ctop%7D "\pmb{\beta}=(\beta_1, \dots, \beta_k)^{\top}"):
  Vetor com ![k](https://latex.codecogs.com/svg.image?k "k")
  coeficientes de regressão, onde
  ![\beta_j](https://latex.codecogs.com/svg.image?%5Cbeta_j "\beta_j"),
  ![j=1, \dots, k](https://latex.codecogs.com/svg.image?j%3D1%2C%20%5Cdots%2C%20k "j=1, \dots, k"),
  indica como a covariável
  ![x_j](https://latex.codecogs.com/svg.image?x_j "x_j") afeta o tempo
  de sobrevivência logarítmico. Um valor
  ![\beta_j](https://latex.codecogs.com/svg.image?%5Cbeta_j "\beta_j")
  positivo sugere que a covariável acelera a falha (diminui o tempo de
  sobrevivência), enquanto um valor negativo implica desaceleração;
- ![\nu_i](https://latex.codecogs.com/svg.image?%5Cnu_i "\nu_i"): Termo
  de erro aleatório associado a cada observação.

Com isso, podemos notar que
<h1 align="center">

![T_i = e^{\mathbf{x}\_i^{\top}\pmb{\beta}}T\_{0i} \\\\\rightarrow \\\\\nu_i = \log(T\_{0i}),](https://latex.codecogs.com/svg.image?T_i%20%3D%20e%5E%7B%5Cmathbf%7Bx%7D_i%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7DT_%7B0i%7D%20%5C%20%5C%20%5Crightarrow%20%5C%20%5C%20%5Cnu_i%20%3D%20%5Clog%28T_%7B0i%7D%29%2C "T_i = e^{\mathbf{x}_i^{\top}\pmb{\beta}}T_{0i} \ \ \rightarrow \ \ \nu_i = \log(T_{0i}),")

</h1>

para o tempo de vida não moderado
![T_0](https://latex.codecogs.com/svg.image?T_0 "T_0") distribuído
independentemente de covariáveis (*baseline*), exceto pela sua própria
igualdade. Diferentes distribuições de
![\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu") implicam
distribuições diferentes de
![T_0](https://latex.codecogs.com/svg.image?T_0 "T_0") (GEORGE; SEALS;
ABAN, 2014).

Com essas definições, pode-se denotar o *AFT model* em termos de
funções:

**Survival function**

A função de sobrevivência de
![T\|\mathbf{x}](https://latex.codecogs.com/svg.image?T%7C%5Cmathbf%7Bx%7D "T|\mathbf{x}")
pode ser descrita pela função de sobrevivência de
![T_0](https://latex.codecogs.com/svg.image?T_0 "T_0"):
<h1 align="center">

![\begin{align\*} 
\mathcal{S}(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = P\_{\pmb{\vartheta}, \pmb{\beta}}(T \> t \\\|\mathbf{x}) \\
& = P\_{\pmb{\vartheta}}\left(e^{\mathbf{x}^{\top}\pmb{\beta}}e^{\nu} \> t \right) \\
& = P\_{\pmb{\vartheta}}\left(T_0 \> te^{-\mathbf{x}^{\top}\pmb{\beta}} \right) \\
& = \mathcal{S}\_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right),
\end{align\*}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%2A%7D%20%0A%5Cmathcal%7BS%7D%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%26%20%3D%20P_%7B%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%7D%28T%20%3E%20t%20%5C%20%7C%5Cmathbf%7Bx%7D%29%20%5C%5C%0A%26%20%3D%20P_%7B%5Cpmb%7B%5Cvartheta%7D%7D%5Cleft%28e%5E%7B%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7De%5E%7B%5Cnu%7D%20%3E%20t%20%5Cright%29%20%5C%5C%0A%26%20%3D%20P_%7B%5Cpmb%7B%5Cvartheta%7D%7D%5Cleft%28T_0%20%3E%20te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%20%5Cright%29%20%5C%5C%0A%26%20%3D%20%5Cmathcal%7BS%7D_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%2C%0A%5Cend%7Balign%2A%7D "\begin{align*} 
\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = P_{\pmb{\vartheta}, \pmb{\beta}}(T > t \ |\mathbf{x}) \\
& = P_{\pmb{\vartheta}}\left(e^{\mathbf{x}^{\top}\pmb{\beta}}e^{\nu} > t \right) \\
& = P_{\pmb{\vartheta}}\left(T_0 > te^{-\mathbf{x}^{\top}\pmb{\beta}} \right) \\
& = \mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right),
\end{align*}")

</h1>

sendo
![\pmb{\vartheta}](https://latex.codecogs.com/svg.image?%5Cpmb%7B%5Cvartheta%7D "\pmb{\vartheta}")
um vetor de parâmetros associado à distribuição *baseline*. Dessa forma,
é possível expressar as outras funções para o *AFT model*.

**Density function**

<h1 align="center">

![\begin{align\*} 
f(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\frac{\partial}{\partial t}\mathcal{S}(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\frac{\partial}{\partial t}\mathcal{S}\_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right) \\
& = -\frac{\partial}{\partial t}\left\[1-F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right) \right\] \\
& = \frac{\partial}{\partial t}F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right) \\
& = f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}},
\end{align\*}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%2A%7D%20%0Af%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%26%20%3D%20-%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5Cmathcal%7BS%7D%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%5C%5C%0A%26%20%3D%20-%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5Cmathcal%7BS%7D_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%20%5C%5C%0A%26%20%3D%20-%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5Cleft%5B1-F_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%20%5Cright%5D%20%5C%5C%0A%26%20%3D%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7DF_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%20%5C%5C%0A%26%20%3D%20f_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29e%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%2C%0A%5Cend%7Balign%2A%7D "\begin{align*} 
f(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\frac{\partial}{\partial t}\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\frac{\partial}{\partial t}\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = -\frac{\partial}{\partial t}\left[1-F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \right] \\
& = \frac{\partial}{\partial t}F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}},
\end{align*}")

</h1>

com
![F_0(\cdot\|\pmb{\vartheta})](https://latex.codecogs.com/svg.image?F_0%28%5Ccdot%7C%5Cpmb%7B%5Cvartheta%7D%29 "F_0(\cdot|\pmb{\vartheta})")
a função de distribuição *baseline*.

**Hazard function**

<h1 align="center">

![\begin{align\*} 
h(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = \frac{f(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x})}{\mathcal{S}(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x})} \\
& = \frac{f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}}{\mathcal{S}\_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right)} \\
& = h_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}
\end{align\*}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%2A%7D%20%0Ah%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%26%20%3D%20%5Cfrac%7Bf%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%7D%7B%5Cmathcal%7BS%7D%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%7D%20%5C%5C%0A%26%20%3D%20%5Cfrac%7Bf_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29e%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%7D%7B%5Cmathcal%7BS%7D_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%7D%20%5C%5C%0A%26%20%3D%20h_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29e%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%0A%5Cend%7Balign%2A%7D "\begin{align*} 
h(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = \frac{f(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x})}{\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x})} \\
& = \frac{f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}}{\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)} \\
& = h_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}
\end{align*}")

</h1>

**Cumulative hazard function**

<h1 align="center">

![\begin{align\*} 
H(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\log\mathcal{S}(t\| \\\pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\log\mathcal{S}\_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right) \\
& = H_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big\|\pmb{\vartheta} \right).
\end{align\*}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%2A%7D%20%0AH%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%26%20%3D%20-%5Clog%5Cmathcal%7BS%7D%28t%7C%20%5C%20%5Cpmb%7B%5Cvartheta%7D%2C%20%5Cpmb%7B%5Cbeta%7D%2C%20%5Cmathbf%7Bx%7D%29%20%5C%5C%0A%26%20%3D%20-%5Clog%5Cmathcal%7BS%7D_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29%20%5C%5C%0A%26%20%3D%20H_0%5Cleft%28te%5E%7B-%5Cmathbf%7Bx%7D%5E%7B%5Ctop%7D%5Cpmb%7B%5Cbeta%7D%7D%5CBig%7C%5Cpmb%7B%5Cvartheta%7D%20%5Cright%29.%0A%5Cend%7Balign%2A%7D "\begin{align*} 
H(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\log\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\log\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = H_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right).
\end{align*}")

</h1>

Os modelos AFT são uma ferramenta valiosa para analisar dados quando há
covariáveis omitidas ou quando a distribuição de probabilidade
subjacente é incerta (BURZYKOWSKI, 2022).

**References:**

- BURZYKOWSKI, Tomasz. Semi‐parametric accelerated failure‐time model: A
  useful alternative to the proportional‐hazards model in cancer
  clinical trials. **Pharmaceutical Statistics**, v. 21, n. 2,
  p. 292-308, 2022.
- GEORGE, Brandon; SEALS, Samantha; ABAN, Inmaculada. Survival analysis
  and regression models. **Journal of nuclear cardiology**, v. 21, n. 4,
  p. 686-694, 2014.
- KALBFLEISCH, John D.; PRENTICE, Ross L. **The statistical analysis of
  failure time data**. John Wiley & Sons, 2002.
