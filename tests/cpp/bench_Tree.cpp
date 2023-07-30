/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 17:45:08 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 21:53:50 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "geoslib_define.h"

#include "Tree/ball.hpp"
#include "Tree/json.hpp"

void	process_object(t_data *data, json_value* value)
{
	int	length;
	int	x;
	int	y;
	int z;
	char *name;
	json_value *temp;

	if (value == NULL)
		return;
	length = value->u.object.length;
	for (x = 0; x < length; x++)
	{
		name = strdup(value->u.object.values[x].name);
		if (!strcmp(name, "n_neighbors"))
			data->n_neighbors = value->u.object.values[x].value->u.integer;
		else if (!strcmp(name, "leaf_size"))
			data->leaf_size = value->u.object.values[x].value->u.integer;
		else if (!strcmp(name, "_tree") || !strcmp(name, "_test"))
		{
			if (!strcmp(name, "_tree"))
				temp = value->u.object.values[x].value->u.object.values[1].value;
			else
				temp = value->u.object.values[x].value;
			data->n_samples = temp->u.array.length;
			data->data = (double**)malloc(sizeof(double*) * data->n_samples);
			for (y =0; y < data->n_samples; y++)
			{
				data->n_features = temp->u.array.values[y]->u.array.length;
				data->data[y] = (double*)malloc(sizeof(double) * data->n_features);
				for (z = 0; z < data->n_features; z++)
					data->data[y][z] = temp->u.array.values[y]->u.array.values[z]->u.dbl;
			}
		}
		free(name);
	}
}

void	process_value(t_data *data, json_value* value)
{
	int j;
	(void)j;
	if (value == NULL)
		return;
	switch (value->type)
	{
		case json_none:
			break;
		case json_object:
			process_object(data, value);
			break;
		case json_array:
			break;
		case json_integer:
			break;
		case json_double:
			break;
		case json_string:
			break;
		case json_boolean:
			break;
		case json_null:
			break;
	}
}

t_data	*read_input(const char *filename)
{
	FILE 		*fp;
	struct stat filestatus;
	int			file_size;
	char		*file_contents;
	json_char	*json;
	json_value	*value;
	t_data		*data;

	if (stat(filename, &filestatus) != 0)
	{
		fprintf(stderr, "File %s not found\n", filename);
		exit(-1);
	}
	file_size = filestatus.st_size;
	file_contents = (char*)malloc(sizeof(char) * filestatus.st_size);
	if (file_contents == NULL)
	{
		fprintf(stderr, "Memory error: unable to allocate %d bytes\n", file_size);
		exit(-1);
	}

	fp = fopen(filename, "rt");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open %s\n", filename);
		free(file_contents);
		exit(-1);
	}
	if (fread(file_contents, file_size, 1, fp) != 1 ) {
		fprintf(stderr, "Unable to read content of %s\n", filename);
		fclose(fp);
		free(file_contents);
		exit(-1);
	}
	fclose(fp);

	json = (json_char*)file_contents;
	value = json_parse(json, file_size);
	if (value == NULL)
	{
		fprintf(stderr, "Unable to parse data\n");
		free(file_contents);
		exit(-11);
	}
	data = (t_data*)malloc(sizeof(t_data));
	process_value(data, value);
	json_value_free(value);
	free(file_contents);
	return (data);
}

t_data	*read_test_data(const char *filename)
{
	FILE		*fp;
	struct stat	filestatus;
	int			file_size;
	char		*file_contents;
	json_char	*json;
	json_value	*value;
	t_data		*data;

	if (stat(filename, &filestatus) != 0)
	{
		fprintf(stderr, "File %s not found\n", filename);
		exit(-1);
	}
	file_size = filestatus.st_size + 12;
	file_contents = (char*)malloc(sizeof(char) * file_size);
	memcpy(file_contents, "{\"_test\": ", 10);
	if (file_contents == NULL)
	{
		fprintf(stderr, "Memory error: unable to allocate %d bytes\n", file_size);
		exit(-1);
	}

	fp = fopen(filename, "rt");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open %s\n", filename);
		free(file_contents);
		exit(-1);
	}
	if (fread(file_contents + 10, file_size - 12, 1, fp) != 1 ) {
		fprintf(stderr, "Unable to read content of %s\n", filename);
		fclose(fp);
		free(file_contents);
		exit(-1);
	}
	fclose(fp);
	file_contents[file_size - 2] = '}';
	file_contents[file_size - 1] = '\0';
	
	json = (json_char*)file_contents;
	value = json_parse(json, file_size);
	if (value == NULL)
	{
		fprintf(stderr, "Unable to parse data\n");
		free(file_contents);
		exit(-11);
	}
	data = (t_data*)malloc(sizeof(t_data));
	process_value(data, value);
	json_value_free(value);
	free(file_contents);
	return (data);
}

void	write_output(t_knn knn, const char *filename)
{
	FILE	*fp;
	int		i, j;

	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Unable to create %s\n", filename);
		exit(-1);
	}
	
	fprintf(fp, "[");
	for (i = 0; i < knn.n_samples; i++)
	{
		fprintf(fp, "{\"knn_idx\": [");
		for (j = 0; j < knn.n_neighbors; j++)
		{
			fprintf(fp, "%d", knn.indices[i][j]);
			if (j != knn.n_neighbors - 1)
				fprintf(fp, ", ");
		}
		fprintf(fp, "], \"knn_dist\": [");
		for (j = 0; j < knn.n_neighbors; j++)
		{
			fprintf(fp, "%.12lf", knn.distances[i][j]);
			if (j != knn.n_neighbors - 1)
				fprintf(fp, ", ");
		}
		fprintf(fp, "]}");
		if (i != knn.n_samples - 1)
			fprintf(fp, ", ");
	}
	fprintf(fp, "]");
	fclose(fp);
}

t_btree	*btree_init_wrapper(t_data *input_data)
{
	return (btree_init(input_data->data, input_data->n_samples, input_data->n_features, input_data->leaf_size));
}

t_knn	btree_query_wrapper(t_btree *tree, t_data *input_data, t_data *test_data)
{
	return (btree_query(tree, test_data->data, test_data->n_samples, test_data->n_features, input_data->n_neighbors)); 
}

void	free_data(t_data *input_data, t_data *test_data)
{
	free_2d_double(input_data->data, input_data->n_samples);
	free(input_data);
	free_2d_double(test_data->data, test_data->n_samples);
	free(test_data);
}

int	main(int argc, char **argv)
{
	t_data	*input_data;
	t_data	*test_data;
	t_btree	*tree;
	t_knn	knn;
	clock_t	start, end;

	String file_param = "/home/drenard/Téléchargements/BallTree-master/datasets/knn_params_6_pts_2_dims.json";
  String file_data = "/home/drenard/Téléchargements/BallTree-master/datasets/test_data_6_pts_2_dims.json";
  String file_results = "/home/drenard/Téléchargements/BallTree-master/datasets/my2_knn_results_6_pts_2_dims.json";

	start = clock();
	input_data = read_input(file_param.c_str());
	end = clock();
	printf("Reading input param took %lf seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

	start = clock();
	test_data = read_test_data(file_data.c_str());
	end = clock();
	printf("Reading test data took %lf seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

	start = clock();
	tree = btree_init_wrapper(input_data);
	end = clock();
	printf("Building ball tree of %d points of %d dimension with leaf size %d took %lf seconds\n",
		input_data->n_samples, input_data->n_features, input_data->leaf_size,
	(double)(end - start) / CLOCKS_PER_SEC);

	start = clock();
	knn = btree_query_wrapper(tree, input_data, test_data);
	end = clock();
	printf("Querying the %d nearest_neighbors of %d points of %d dimension took %lf seconds\n",
		input_data->n_neighbors, test_data->n_samples, test_data->n_features, (double)(end - start) / CLOCKS_PER_SEC);

	start = clock();
	write_output(knn, file_results.c_str());
	end = clock();
	printf("Writing result data took %lf second\n", (double)(end - start) / CLOCKS_PER_SEC);

	//free stuff
	free_tree(tree);
	free_knn(knn, test_data->n_samples);
	free_data(input_data, test_data);
	return (0);
}
