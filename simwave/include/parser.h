#ifndef _PARSER_H_
#define _PARSER_H_
///
/// \copyright Copyright 2015 Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// \b simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// \author Issam Said
/// \file parser.h
/// \version $Id: parser.h 2815 2015-10-27 09:42:35Z isaid $
/// \brief A customized user options parser.
///
/// Parses user options from command line or from a text file.
///

#include <stdlib.h>



/// \brief Describes the type of a user option
enum _option_kind {
    BOOL     = 0,
    INT         ,
    FLOAT       ,
    STRING      ,
    VECT_INT    ,
    VECT_FLOAT  ,
    VECT_STRING
};

/// \brief Defines a boolean type
typedef char bool_t;

/// \brief Typedef \ref _option_kind for convenience
typedef enum _option_kind option_kind;

/// \brief Typedef \ref _option for convenience
typedef struct _option option;

/// \brief Typedef \ref _parser for convenience
typedef struct _parser parser;
/// \brief A struct that describes a single user option
struct _option {
    option_kind kind;
    char    short_name;
    char   *long_name;
    char   *default_value;
    char   *value;
    char   *message;
    option *next;
};

/// \brief A linked list that contains all the user options
struct _parser {
    char exe_name[128];
    char* description;
    option* options;
};

/// \brief Creates a user option
/// \param kind the type of the option
/// \param short_name is a single character to reference the option
/// \param long_name a long name that characterizes the option
/// \param default_value is the default value assigned to the option once created
/// \param message a help message that describes the option
/// \return a pointer to the option
option* option_create(const option_kind kind,
                      const char short_name,
                      const char *long_name,
                      const char *default_value,
                      const char *message);

/// \brief Deletes a user option
/// \param o pointer to the option to be deleted
void option_delete(option *o);

/// \brief Prints informations about a user option
/// \param o pointer to the option to be printed
void option_print(option *o);

/// \brief Checks if a user option is of \ref BOOL kind
/// \param o pointer to the option to be checked
bool_t option_boolean(option *o);

/// \brief Deletes a user option
/// \param o pointer to the option to be deleted
void option_delete(option *o);

/// \brief Prints informations about a user option
/// \param o pointer to the option to be printed
void option_print(option *o);

/// \brief Checks if a user option is of a \ref VECT_INT,
/// of a \ref VECT_STRING or of a \ref VECT_FLOAT kind
/// \param o pointer to the option to be checked
bool_t option_vect(option *o);

/// \brief Creates a parser
/// \return a pointer to the created parser
parser* parser_create();

/// \brief Deletes a parser
/// \param p pointer to the parser to be deleted
void parser_delete(parser* p);

/// \brief Prints all available options to the user
/// \param p pointer to the parser to be printed
void parser_usage(parser* p);

/// \brief Finds an option by its long name
/// \param p pointer to the parser that may contain the option
/// \param name the long name of the option
/// \return a pointer to he option if found or NULL
option* parser_find_option_by_long_name(parser* p, const char* name);

/// \brief Finds an option by its short name
/// \param p pointer to the parser that may contain the option
/// \param name the short name of the option
/// \return a pointer to he option if found or NULL
option* parser_find_option_by_short_name(parser* p, char name);

/// \brief Parses the command line entries
/// \param p pointer to a parser that will contain the parsed options
/// \param argc the number of command line entries
/// \param argv an array that contains the command line entries
/// \return the index of the last parsed entry
int  parser_parse(parser* p, int argc, char* argv[]);

/// \brief Looks for options in a text file.
/// The format of the text entries should be key=value
/// \param p pointer to a parser that will contain the parsed options
/// \param file_name the name of the text file to be parsed
/// \return the index of the last parsed entry
int  parser_parse_from_file(parser* p, char* file_name);

/// \brief Creates and adds a new option to an existing parser
/// \param p a pointer to the parser
/// \param kind the type of the option
/// \param short_name is a single character to reference the option
/// \param long_name a long name that characterizes the option
/// \param default_value is the default value assigned to the option once created
/// \param message a help message that describes the option
void parser_put(parser* p,
                const option_kind kind,
                const char short_name,
                const char *long_name,
                const char *default_value,
                const char *message);

/// \brief Checks if the user has asked help
/// \param p a pointer to the parser
bool_t  parser_help(parser* p);

/// \brief Detects the value of a \ref BOOL option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \return the value of the option
bool_t  parser_get_bool(parser* p, const char* expression);

/// \brief Detects the value of a \ref INT option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \return the value of the option
int   parser_get_int(parser* p, const char* expression);

/// \brief Detects the value of a \ref FLOAT option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \return the value of the option
float parser_get_float(parser* p, const char* expression);

/// \brief Detects the value of a \ref STRING option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \return the value of the option
char* parser_get_string(parser* p, const char* expression);

/// \brief Detects the values of a \ref VECT_INT option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \param v a pointer to an array that will contains the values of the option
/// \param n the number of the values of the option
void  parser_get_vect_int(parser* p, const char* expression, int* v, unsigned int n);

/// \brief Detects the values of a \ref VECT_INT option (size_t)
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \param v a pointer to an array that will contains the values of the option
/// \param n the number of the values of the option
void  parser_get_vect_size_t(parser* p, const char* expression, size_t* v, unsigned int n);

/// \brief Detects the values of a \ref VECT_FLOAT option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \param vf a pointer to an array that will contains the values of the option
/// \param n the number of the values of the option
void  parser_get_vect_float(parser* p, const char* expression, float* vf, unsigned int n);

/// \brief Detects the values of a \ref VECT_STRING option
/// \param p a pointer to the parser
/// \param expression a string that contains the long name or
///  the short name of the option
/// \param vs a pointer to an array that will contains the values of the option
/// \param n the number of the values of the option
void  parser_get_vect_string(parser* p, const char* expression, char** vs, unsigned int n);

/// \brief Sets the possible options for \b simwave
///
/// run <b> bin/simwave -h </b> or <b> bin/simwave --help </b> to see thie possible options
#define PARSER_BOOTSTRAP(p)                                                                  \
parser_put(p, BOOL,     'v', "verbose",    "false",         "enable verbose mode");          \
parser_put(p, BOOL,     'm', "ms",         "false",         "use milliseconds for timing");  \
parser_put(p, BOOL,     'u', "us",         "false",         "use microseconds for timing");  \
parser_put(p, BOOL,     'n', "ns",         "false",         "use nanoseconds  for timing");  \
parser_put(p, VECT_INT,   0, "grid",       "100,100,100",   "domain x dimension");           \
parser_put(p, VECT_INT,   0, "dgrid",      "10,10,10",      "space delta step");             \
parser_put(p, VECT_INT,   0, "source_loc", "-1,-1,-1",      "source location in the grid");  \
parser_put(p, INT,      'i', "iter",      "1000",           "simulation time step number");  \
parser_put(p, FLOAT,      0, "cfl",        "0.8",           "cfl pourcentage");              \
parser_put(p, FLOAT,      0, "fmax",       "25.",           "source max frequency");         \
parser_put(p, FLOAT,      0, "vmin",       "1500.",         "min velocity");                 \
parser_put(p, FLOAT,      0, "vmax",       "4500.",         "max velocity");                 \
parser_put(p, STRING,     0, "output",     "data/wave.bin", "snapshots file path");          \
parser_put(p, BOOL,       0, "nosnap",     "false",         "do not save snapshots");        \
parser_put(p, INT,        0, "nbsnap",     "10",             "snapshot frequency");           \
parser_put(p, VECT_INT,   0, "bloc",       "32,4,2",   "x * y * z");           \
parser_put(p, BOOL,     0, "opt",    "false",         "enable optimal mode");          \
parser_put(p, BOOL,       0, "check",      "false",         "check the gpu results");        \
parser_put(p, FLOAT,    'e', "epsilon",    "1.e-2", "the margin of floating point errors")

/// \brief Adds the PML options to the parser
#define PARSER_ADD_PML_OPTIONS(p) parser_put(p, VECT_INT, 0, "taper", "2,2,2", "pml tapers")

#endif //  _PARSER_H_

