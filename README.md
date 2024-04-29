Accelerated survival analysis
================

<style>
body {
text-align: justify}
</style>

Os modelos de regressão de sobrevivência de classe acelerada são um
grupo de técnicas estatísticas usadas para analisar dados em que a
variável de resultado é o tempo até um evento, também conhecidos como
dados de sobrevivência. Estes modelos diferem dos modelos de
sobrevivência padrão por incorporarem o efeito das covariáveis na
própria escala de tempo, em vez de apenas na taxa de risco. Aqui está
uma análise dos três tipos principais (AFT, AH, AO):

## 1 - Modelo de tempo de falha acelerado (*AFT model*)

O AFT model é uma abordagem paramétrica (ou semiparamétrica) que oferece
uma opção robusta de regressão em análise de sobrevivência. A ideia
central por trás dos modelos AFT é que as covariáveis atuam
multiplicativamente no log do tempo até a ocorrência de algum evento de
interesse, acelerando ou desacelerando o processo de falha. Assim, para
a $i$-ésima observação, um AFT model é caracterizado por sugerir que $$
\log(T_i) = \mathbf{x}_i^{\top}\pmb{\beta} + \nu_i, \ \ \ \ i=1, \dots, n,
$$ em que

- $T_i$: $i$-ésima observação da variável aleatória que define o tempo
  de falha ($T>0$);
- $\mathbf{x}=(x_1, \dots, x_k)^{\top}$: Vetor de covariáveis com
  dimensão $k$ para a $i$-ésima observação;
- $\pmb{\beta}=(\beta_1, \dots, \beta_k)^{\top}$: Vetor com $k$
  coeficientes de regressão, onde $\beta_j$, $j=1, \dots, k$, indica
  como a covariável $x_j$ afeta o tempo de sobrevivência logarítmico. Um
  valor $\beta_j$ positivo sugere que a covariável acelera a falha
  (diminui o tempo de sobrevivência), enquanto um valor negativo implica
  desaceleração;
- $\nu_i$: Termo de erro aleatório associado a cada observação.

Com isso, podemos notar que

$$
T_i = e^{\mathbf{x}_i^{\top}\pmb{\beta}}e^{\nu_i} \ \ \rightarrow \ \ T_{0i} = e^{\nu_i} = T_ie^{-\mathbf{x}_i^{\top}\pmb{\beta}},
$$ para o tempo de vida não moderado $T_0$ distribuído independentemente
de covariáveis (baseline), exceto pela sua própria igualdade. Diferentes
distribuições de $\nu$ implicam distribuições diferentes de $T_0$.

Com essas definições, pode-se denotar o AFT model em termos de funções:

**Survival function**

A função de sobrevivência de $T|\mathbf{x}$ pode ser descrita pela
função de sobrevivência de $T_0$: $$
\begin{align*} 
\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = P_{\pmb{\vartheta}, \pmb{\beta}}(T > t \ |\mathbf{x}) \\
& = P_{\pmb{\vartheta}}\left(e^{\mathbf{x}^{\top}\pmb{\beta}}e^{\nu} > t \right) \\
& = P_{\pmb{\vartheta}}\left(T_0 > te^{-\mathbf{x}^{\top}\pmb{\beta}} \right) \\
& = \mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right),
\end{align*}
$$ sendo $\pmb{\vartheta}$ um vetor de parâmetros associado à
distribuição baseline. Dessa forma, é possível expressar as outras
funções para o AFT model.

**Density function** $$
\begin{align*} 
f(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\frac{\partial}{\partial t}\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\frac{\partial}{\partial t}\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = -\frac{\partial}{\partial t}\left[1-F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \right] \\
& = \frac{\partial}{\partial t}F_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}},
\end{align*}
$$ com $F_0(\cdot|\pmb{\vartheta})$ a função de distribuição baseline.

**Hazard function** $$
\begin{align*} 
h(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = \frac{f(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x})}{\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x})} \\
& = \frac{f_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}}{\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)} \\
& = h_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right)e^{-\mathbf{x}^{\top}\pmb{\beta}}
\end{align*}
$$

**Cumulative hazard function** $$
\begin{align*} 
H(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) & = -\log\mathcal{S}(t| \ \pmb{\vartheta}, \pmb{\beta}, \mathbf{x}) \\
& = -\log\mathcal{S}_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right) \\
& = H_0\left(te^{-\mathbf{x}^{\top}\pmb{\beta}}\Big|\pmb{\vartheta} \right).
\end{align*}
$$

As estimativas dos parâmetros de regressão dos modelos AFT são robustas
para covariáveis omitidas. Eles também são menos afetados pela escolha
da distribuição de probabilidade. Os resultados dos modelos AFT são
facilmente interpretados. Por exemplo, os resultados de um ensaio
clínico com mortalidade como ponto final poderiam ser interpretados como
um certo aumento percentual na expectativa de vida futura com o novo
tratamento em comparação com o controle. Assim, um paciente poderia ser
informado de que se esperaria que ele vivesse (digamos) 15% mais se
fizesse o novo tratamento. As taxas de risco podem ser mais difíceis de
explicar em termos leigos.

**References:**

- KALBFLEISCH, John D.; PRENTICE, Ross L. **The statistical analysis of
  failure time data**. *John Wiley & Sons*, 2011.
- For a more technical treatment, you can refer to research articles on
  AFT models with specific distributions, such as:
  - On estimation for accelerated failure time models with small or rare
    event survival data:
    <https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-022-01638-1>

I hope this explanation provides a good foundation for understanding the
mathematical concepts and applications of AFT models in survival
analysis.
