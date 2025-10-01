/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <leonardo.robol@unipi.it>
*/

/**
* @file
* @brief Generic input-output functions inside MPSolve.
*/

#ifndef MPS_INPUT_OUTPUT_H_
#define MPS_INPUT_OUTPUT_H_

MPS_BEGIN_DECLS

    void mps_readroots(mps_context* ctx);
void mps_countroots(mps_context* ctx);
void mps_outroot(mps_context* ctx, int i, int num);
void mps_output(mps_context* ctx);
void mps_copy_roots(mps_context* ctx);
void mps_dump_status(mps_context* ctx, FILE* outstr);
void mps_dump(mps_context* ctx);
void mps_dump_cluster_structure(mps_context* ctx, FILE* outstr);
mps_boolean mps_is_a_tty(FILE* stream);
void mps_warn(mps_context* st, char* format, ...);
void mps_error(mps_context* st, const char* format, ...);
void mps_print_errors(mps_context* ctx);

void mps_skip_comments(FILE* input_stream);
void mps_raise_parsing_error(mps_context* ctx, mps_input_buffer* buffer,
    const char* token,
    const char* message, ...);
mps_input_option mps_parse_option_line(mps_context* ctx, char* line, size_t length);

mps_polynomial* mps_monomial_poly_read_from_stream_v2(mps_context* ctx, mps_input_buffer* buffer);

mps_polynomial* mps_monomial_yacc_parser(mps_context* ctx, mps_abstract_input_stream* stream);


MPS_END_DECLS

#endif /* MPS_INPUT_OUTPUT_H_ */
