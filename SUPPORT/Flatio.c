#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define TRUE 1
#define FALSE 0
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct data_format {
	int length;
	char *nuc;
	int offset;
	char name[64];
	char type;
};

char *Realloc(char *block, int size);
char *Calloc(int count, int size);
int ErrorOut(int code, char *string);
int Errorout(char *string);
int ReadFlat(FILE *file, struct data_format align[], int maxseqs);
int WriteData(FILE *file, struct data_format data[], int count);

int ReadFlat(FILE *file, struct data_format align[], int maxseqs)
{
	int j, len = 0, count = -1, offset;
	unsigned maxlen = 1024;
	char cinline[1025];
	extern char *Calloc(), *Realloc();

	if (file == NULL) Errorout("Cannot open data file");

	for (; fgets(cinline, 1024, file) != NULL;) {
		cinline[strlen(cinline) - 1] = '\0';
		switch (cinline[0]) {
			case '>':
			case '#':
			case '%':
			case '"':
			case '@':
				offset = 0;
				for (j = 0; j < strlen(cinline); j++) {
					if (cinline[j] == '(') {
						sscanf(
						    (char *)(cinline + j + 1),
						    "%d", &offset);
						cinline[j] = '\0';
					}
				}

				if (count != -1) {
					align[count].length = len;
					align[count].nuc[len] = '\0';
					maxlen = len;
				}

				count++;
				if (count > maxseqs)
					Errorout(
					    "Sorry, alignment is too large");

				align[count].nuc = Calloc(maxlen, sizeof(char));
				align[count].type = cinline[0];
				align[count].offset = offset;
				if (align[count].nuc == NULL)
					Errorout("Calloc problem");

				sscanf((char *)(cinline + 1), "%s",
				       align[count].name);
				len = 0;
				break;
			default:
				if (len + strlen(cinline) > maxlen) {
					maxlen = (maxlen + strlen(cinline)) * 2;
					align[count].nuc =
					    Realloc(align[count].nuc, maxlen);
				}
				for (j = 0; j < strlen(cinline); j++)
					align[count].nuc[j + len] = cinline[j];
				len += strlen(cinline);
				break;
		}
	}
	if (count == -1) exit(1);

	align[count].length = len;
	align[count].nuc[len] = '\0';
	return (++count);
}

int Errorout(char *string)
{
	fprintf(stderr, "%s\n", string);
	exit(1);
}

int WriteData(FILE *file, struct data_format data[], int count)
{
	int i, j;
	for (j = 0; j < count; j++) {
		if (data[j].offset)
			fprintf(file, "\n%c%s(%d)", data[j].type, data[j].name,
				data[j].offset);
		else
			fprintf(file, "\n%c%s", data[j].type, data[j].name);

		for (i = 0; i < data[j].length; i++) {
			if (i % 60 == 0) fputc('\n', file);
			fputc(data[j].nuc[i], file);
		}
	}
	return 0;
}

int ErrorOut(int code, char *string)
{
	if (code == 0) {
		fprintf(stderr, "Error:%s\n", string);
		exit(1);
	}
	return 0;
}

char *Calloc(int count, int size)
{
	char *temp;

	temp = (char *)calloc(count, size);
	if (temp == NULL) {
		fprintf(stdout, "Error in Calloc\n");
		exit(-1);
	}
	else
		return (temp);
}

char *Realloc(char *block, int size)
{
	char *temp;
	temp = (char *)realloc(block, size);
	if (temp == NULL) {
		fprintf(stdout, "Error in Calloc\n");
		exit(-1);
	}
	else
		return (temp);
}

