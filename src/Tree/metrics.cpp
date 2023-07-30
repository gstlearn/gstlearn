/* MIT License
 *
 * Copyright (C) 2012, 2013, 2014 James McLaughlin et al.  All rights reserved.
 * https://github.com/udp/json-parser
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metrics.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:32 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 16:52:59 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "Tree/ball.hpp"

double	manhattan_dist(double *x1, double *x2, int size)
{
	double	d = 0;

	for (int i = 0; i < size; i++)
		d += fabs(x1[i] - x2[i]);
	return (d);
}

double	min_dist(t_btree *tree, int i_node, double *pt)
{
	double	dist_pt;

	dist_pt = manhattan_dist(pt, tree->node_bounds[0][i_node], tree->n_features);
	return (fmax(0.0, dist_pt - tree->node_data[i_node].radius));
}
