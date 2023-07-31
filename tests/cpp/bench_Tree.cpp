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

#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/AStringable.hpp"
#include "Tree/Ball.hpp"
#include "Tree/ball.h"
#include "Tree/json.h"

static void process_object(t_data *data, json_value* value)
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

static void process_value(t_data *data, json_value* value)
{
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

static t_data *read_input_data(const char *filename)
{
	FILE 		    *fp;
	struct stat  filestatus;
	int			     file_size;
	char		    *file_contents;
	json_char 	*json;
	json_value	*value;
	t_data		  *data;

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

static t_data *read_test_data(const char *filename)
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

static void write_output(t_knn knn, const char *filename)
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

static void free_data(t_data *input_data, t_data *test_data)
{
	free_2d_double(input_data->data, input_data->n_samples);
	free(input_data);
	free_2d_double(test_data->data, test_data->n_samples);
	free(test_data);
}

int main(int argc, char *argv[])
{
	t_data	*input_data = nullptr;
	t_data	*test_data = nullptr;
	t_knn	knn;
	t_knn knn1;
  Timer timer;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

	String file_param = "/home/drenard/Téléchargements/BallTree-master/datasets/knn_params_10k_pts_512_dims.json";
	String file_data  = "/home/drenard/Téléchargements/BallTree-master/datasets/test_data_10k_pts_512_dims.json";
	String file_results = "/home/drenard/Téléchargements/BallTree-master/datasets/my_knn_results_10k_pts_512_dims.json";

  timer.reset();
	input_data = read_input_data(file_param.c_str());
  timer.displayIntervalMilliseconds("Reading Input Parameters and Data", 0);

  timer.reset();
  Ball ball(input_data->data, input_data->n_samples, input_data->n_features, input_data->leaf_size, 0);
  timer.displayIntervalMilliseconds("Instantiating Ball class", 0);

  timer.reset();
  test_data = read_test_data(file_data.c_str());
  timer.displayIntervalMilliseconds("Reading Test Data", 0);

  timer.reset();
  (void) ball.build();
  message("Building ball tree of %d points of dimension %d with leaf size %d\n",
          input_data->n_samples, input_data->n_features, input_data->leaf_size);
	timer.displayIntervalMilliseconds("Building ball tree",0);

  message("\nQuerying the %d nearest_neighbors for %d points of dimension %d\n",
          input_data->n_neighbors, test_data->n_samples, test_data->n_features);

  timer.reset();
	knn = ball.query(test_data->data, test_data->n_samples, test_data->n_features, input_data->n_neighbors);
	display(knn, 3); // Display of the information for the first three tested samples
  timer.displayIntervalMilliseconds("Querying the Ball Tree",0);

	for (int is = 0; is < 3; is++)
	{
	  timer.reset();
	  knn1 = ball.queryOne(test_data->data[is], test_data->n_features, input_data->n_neighbors);
	  display(knn1);
    timer.displayIntervalMilliseconds("Querying the Ball Tree Per Target",0);
	}

	timer.reset();
	write_output(knn, file_results.c_str());
	timer.displayIntervalMilliseconds("Writing result data", 0);

	//free stuff
	free_knn(knn, test_data->n_samples);
	free_data(input_data, test_data);
	return (0);
}
