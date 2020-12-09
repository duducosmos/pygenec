#!/usr/bin/env python3.6
# -*- Coding: UTF-8 -*-
"""
Evolução.

Programa sob licença GNU V.3.
Desenvolvido por: E. S. Pereira.
Versão 0.0.1.
"""

from progressbar import progressbar
from numpy import array
from pathos.helpers import cpu_count, shutdown
from pathos.multiprocessing import ProcessPool

from .evolucao import Evolucao
from .populacao import Populacao


class AlgoritmoGenetico:
    def __init__(self, selecao, cruzamento, mutacao, bits, cromossomos,  **kwargs):
        self._bits = bits
        self._cromossomos = cromossomos

        self._parametros = {'tamanho_populacao': 100,
                            'nsele': 20,
                            'elitista': True,
                            'pcruz': 0.6,
                            'epidemia': 50,
                            'total_geracoes': 100,
                            'pmut': 0.05}

        for key in kwargs:
            if key not in self._parametros:
                raise NameError(f"Parametro invalido {key}")
            self._parametros[key] = kwargs[key]

        genes = self._bits * self._cromossomos 


        self._populacao = Populacao(self.avaliacao,
                                    genes, 
                                    self._parametros['tamanho_populacao'])
        self._cruzamento = cruzamento(self._parametros['tamanho_populacao'])

        self._mutacao = mutacao(pmut=self._parametros['pmut'])

        self._selecao = selecao(self._populacao)
        self._evolucao = Evolucao(self._populacao, self._selecao, self._cruzamento, self._mutacao)

        self._evolucao.nsele = self._parametros['nsele']
        self._evolucao.pcruz = self._parametros['pcruz']
        self._evolucao.manter_melhor = self._parametros['elitista']
        self._evolucao.epidemia = self._parametros['epidemia']

    def funcao_objetivo(self, individuo):
        raise NotImplementedError("Deve ser implementada")

    @property
    def parametros(self):
      return self._parametros

    def avaliacao(self, populacao):

        n = len(populacao)

        def steps(k):
            individuo = populacao[k, :]
            obj = self.funcao_objetivo(individuo)
            return obj

        ncpu = cpu_count()
        pool = ProcessPool(nodes=ncpu)
        pesos = array(pool.map(steps, range(n)))
        pool.close()
        pool.join()
        pool.clear()
        shutdown()

        return pesos

    def __call__(self):

        for _ in progressbar(range(self._parametros['total_geracoes'])):
            vmin, vmax = self._evolucao.evoluir()
            yield self._populacao.populacao[-1, :]