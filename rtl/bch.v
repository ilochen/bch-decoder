module bch(
	input clk,
	input rstn,
	input mode,
	input [1:0] code,
	input set,
	input [63:0] idata,
	output ready,
	output finish,
	output [9:0] odata
);

reg [2:0] state_w, state_r;
localparam IDLE = 3'b000,
			SYNDROME_GEN = 3'b001,
			ERROR_POLY_DISCREPENCY = 3'b010,
			ERROR_POLY_UPDATE = 3'b011,
			ERROR_LOCATE = 3'b100,
			CORRELATION = 3'b101,
			OUTPUT = 3'b110;

reg [1:0] code_w, code_r;
reg mode_w, mode_r;
reg output_mode_w, output_mode_r;
reg ready_w, ready_r;
reg finish_w, finish_r;
reg [7:0] syndrome_data_in_w, syndrome_data_in_r;
reg [7:0] syndrome_data_in_piped_w, syndrome_data_in_piped_r;
reg [7:0] syndrome_data_flipped_w [0:2], syndrome_data_flipped_r [0:2];
wire clear;
assign clear = (state_r != SYNDROME_GEN) ? 1'b1 : 1'b0;


reg [7:0] counter_w, counter_r;
reg [6:0] llr_w [0:1023], llr_r [0:1023]; // max 128 symbols

reg [8:0] input_done_w, input_done_r;

wire [9:0] syndrome_out [0:3][0:7];
reg [9:0] syndrome_w [0:3][0:7], syndrome_r [0:3][0:7]; // max 8 syndromes for BCH(1023,983)


reg [9:0] sigma_w [0:3][0:4], sigma_r [0:3][0:4]; // max degree error polynomial for BCH(1023,983) = 4
reg [9:0] current_syndrome_w [0:3][0:4], current_syndrome_r [0:3][0:4];
reg [9:0] delta_w [0:3], delta_r [0:3]; // discrepancy
reg [9:0] gamma_w [0:3], gamma_r [0:3]; // scaling factor
reg [3:0] L_w [0:3], L_r [0:3]; // degree of error locator polynomial
reg [3:0] m_w [0:3], m_r [0:3]; // iteration counter
reg [9:0] Bx_w [0:3][0:4], Bx_r [0:3][0:4]; // previous error locator polynomial

reg [9:0] gf_mult_in1 [0:3][0:9];
reg [9:0] gf_mult_in2 [0:3][0:9];
wire [9:0] gf_mult_out [0:3][0:9];

wire [9:0] chien_shifted_value [0:3][0:7][0:4];
reg [9:0] error_position_w [0:3][0:3], error_position_r [0:3][0:3]; // max 6 errors for BCH(1023,983) soft decoding
reg [2:0] error_count_w[0:3], error_count_r[0:3];
reg [7:0] error_found_w[0:3], error_found_r[0:3];
reg [2:0] error_count [0:3][0:6];
reg [6:0] partial_sum [0:3][0:7];
reg [7:0] partial_sum_w [0:3][0:3], partial_sum_r [0:3][0:3];
reg init;

reg [9:0] error_output_w, error_output_r;
reg [9:0] corr_w [0:3], corr_r [0:3];
reg [1:0] corr_index_w, corr_index_r;
wire [1:0] corr_index;

always @(*) begin
	error_found_w[0][0] = ((sigma_r[0][0]^ chien_shifted_value[0][7][1] ^ chien_shifted_value[0][7][2] ^ chien_shifted_value[0][7][3] ^ chien_shifted_value[0][7][4]) == 10'd0);
	error_found_w[0][1] = ((sigma_r[0][0]^ chien_shifted_value[0][6][1] ^ chien_shifted_value[0][6][2] ^ chien_shifted_value[0][6][3] ^ chien_shifted_value[0][6][4]) == 10'd0);
	error_found_w[0][2] = ((sigma_r[0][0]^ chien_shifted_value[0][5][1] ^ chien_shifted_value[0][5][2] ^ chien_shifted_value[0][5][3] ^ chien_shifted_value[0][5][4]) == 10'd0);
	error_found_w[0][3] = ((sigma_r[0][0]^ chien_shifted_value[0][4][1] ^ chien_shifted_value[0][4][2] ^ chien_shifted_value[0][4][3] ^ chien_shifted_value[0][4][4]) == 10'd0);
	error_found_w[0][4] = ((sigma_r[0][0]^ chien_shifted_value[0][3][1] ^ chien_shifted_value[0][3][2] ^ chien_shifted_value[0][3][3] ^ chien_shifted_value[0][3][4]) == 10'd0);
	error_found_w[0][5] = ((sigma_r[0][0]^ chien_shifted_value[0][2][1] ^ chien_shifted_value[0][2][2] ^ chien_shifted_value[0][2][3] ^ chien_shifted_value[0][2][4]) == 10'd0);
	error_found_w[0][6] = ((sigma_r[0][0]^ chien_shifted_value[0][1][1] ^ chien_shifted_value[0][1][2] ^ chien_shifted_value[0][1][3] ^ chien_shifted_value[0][1][4]) == 10'd0);
	error_found_w[0][7] = ((sigma_r[0][0]^ chien_shifted_value[0][0][1] ^ chien_shifted_value[0][0][2] ^ chien_shifted_value[0][0][3] ^ chien_shifted_value[0][0][4]) == 10'd0);

	error_found_w[1][0] = ((sigma_r[1][0]^ chien_shifted_value[1][7][1] ^ chien_shifted_value[1][7][2] ^ chien_shifted_value[1][7][3] ^ chien_shifted_value[1][7][4]) == 10'd0);
	error_found_w[1][1] = ((sigma_r[1][0]^ chien_shifted_value[1][6][1] ^ chien_shifted_value[1][6][2] ^ chien_shifted_value[1][6][3] ^ chien_shifted_value[1][6][4]) == 10'd0);
	error_found_w[1][2] = ((sigma_r[1][0]^ chien_shifted_value[1][5][1] ^ chien_shifted_value[1][5][2] ^ chien_shifted_value[1][5][3] ^ chien_shifted_value[1][5][4]) == 10'd0);
	error_found_w[1][3] = ((sigma_r[1][0]^ chien_shifted_value[1][4][1] ^ chien_shifted_value[1][4][2] ^ chien_shifted_value[1][4][3] ^ chien_shifted_value[1][4][4]) == 10'd0);
	error_found_w[1][4] = ((sigma_r[1][0]^ chien_shifted_value[1][3][1] ^ chien_shifted_value[1][3][2] ^ chien_shifted_value[1][3][3] ^ chien_shifted_value[1][3][4]) == 10'd0);
	error_found_w[1][5] = ((sigma_r[1][0]^ chien_shifted_value[1][2][1] ^ chien_shifted_value[1][2][2] ^ chien_shifted_value[1][2][3] ^ chien_shifted_value[1][2][4]) == 10'd0);
	error_found_w[1][6] = ((sigma_r[1][0]^ chien_shifted_value[1][1][1] ^ chien_shifted_value[1][1][2] ^ chien_shifted_value[1][1][3] ^ chien_shifted_value[1][1][4]) == 10'd0);
	error_found_w[1][7] = ((sigma_r[1][0]^ chien_shifted_value[1][0][1] ^ chien_shifted_value[1][0][2] ^ chien_shifted_value[1][0][3] ^ chien_shifted_value[1][0][4]) == 10'd0);

	error_found_w[2][0] = ((sigma_r[2][0]^ chien_shifted_value[2][7][1] ^ chien_shifted_value[2][7][2] ^ chien_shifted_value[2][7][3] ^ chien_shifted_value[2][7][4]) == 10'd0);
	error_found_w[2][1] = ((sigma_r[2][0]^ chien_shifted_value[2][6][1] ^ chien_shifted_value[2][6][2] ^ chien_shifted_value[2][6][3] ^ chien_shifted_value[2][6][4]) == 10'd0);
	error_found_w[2][2] = ((sigma_r[2][0]^ chien_shifted_value[2][5][1] ^ chien_shifted_value[2][5][2] ^ chien_shifted_value[2][5][3] ^ chien_shifted_value[2][5][4]) == 10'd0);
	error_found_w[2][3] = ((sigma_r[2][0]^ chien_shifted_value[2][4][1] ^ chien_shifted_value[2][4][2] ^ chien_shifted_value[2][4][3] ^ chien_shifted_value[2][4][4]) == 10'd0);
	error_found_w[2][4] = ((sigma_r[2][0]^ chien_shifted_value[2][3][1] ^ chien_shifted_value[2][3][2] ^ chien_shifted_value[2][3][3] ^ chien_shifted_value[2][3][4]) == 10'd0);
	error_found_w[2][5] = ((sigma_r[2][0]^ chien_shifted_value[2][2][1] ^ chien_shifted_value[2][2][2] ^ chien_shifted_value[2][2][3] ^ chien_shifted_value[2][2][4]) == 10'd0);
	error_found_w[2][6] = ((sigma_r[2][0]^ chien_shifted_value[2][1][1] ^ chien_shifted_value[2][1][2] ^ chien_shifted_value[2][1][3] ^ chien_shifted_value[2][1][4]) == 10'd0);
	error_found_w[2][7] = ((sigma_r[2][0]^ chien_shifted_value[2][0][1] ^ chien_shifted_value[2][0][2] ^ chien_shifted_value[2][0][3] ^ chien_shifted_value[2][0][4]) == 10'd0);

	error_found_w[3][0] = ((sigma_r[3][0]^ chien_shifted_value[3][7][1] ^ chien_shifted_value[3][7][2] ^ chien_shifted_value[3][7][3] ^ chien_shifted_value[3][7][4]) == 10'd0);
	error_found_w[3][1] = ((sigma_r[3][0]^ chien_shifted_value[3][6][1] ^ chien_shifted_value[3][6][2] ^ chien_shifted_value[3][6][3] ^ chien_shifted_value[3][6][4]) == 10'd0);
	error_found_w[3][2] = ((sigma_r[3][0]^ chien_shifted_value[3][5][1] ^ chien_shifted_value[3][5][2] ^ chien_shifted_value[3][5][3] ^ chien_shifted_value[3][5][4]) == 10'd0);
	error_found_w[3][3] = ((sigma_r[3][0]^ chien_shifted_value[3][4][1] ^ chien_shifted_value[3][4][2] ^ chien_shifted_value[3][4][3] ^ chien_shifted_value[3][4][4]) == 10'd0);
	error_found_w[3][4] = ((sigma_r[3][0]^ chien_shifted_value[3][3][1] ^ chien_shifted_value[3][3][2] ^ chien_shifted_value[3][3][3] ^ chien_shifted_value[3][3][4]) == 10'd0);
	error_found_w[3][5] = ((sigma_r[3][0]^ chien_shifted_value[3][2][1] ^ chien_shifted_value[3][2][2] ^ chien_shifted_value[3][2][3] ^ chien_shifted_value[3][2][4]) == 10'd0);
	error_found_w[3][6] = ((sigma_r[3][0]^ chien_shifted_value[3][1][1] ^ chien_shifted_value[3][1][2] ^ chien_shifted_value[3][1][3] ^ chien_shifted_value[3][1][4]) == 10'd0);
	error_found_w[3][7] = ((sigma_r[3][0]^ chien_shifted_value[3][0][1] ^ chien_shifted_value[3][0][2] ^ chien_shifted_value[3][0][3] ^ chien_shifted_value[3][0][4]) == 10'd0);
end
reg [9:0] odata_r, odata_w;
// output assignments
assign ready = ready_r;
assign finish = finish_r;
assign odata = odata_r;

// for soft decoding
wire [6:0] abs_llr [0:7];
reg [9:0] min_idx0_w, min_idx0_r;
reg [9:0] min_idx1_w, min_idx1_r;
reg [6:0] min_val0_w, min_val0_r;
reg [6:0] min_val1_w, min_val1_r;
wire [9:0] min_idx0, min_idx1;
wire [6:0] min_val0, min_val1;

wire [1:0] fetch_error_pos;
assign fetch_error_pos = error_count_r[corr_index_r] - 3'd1;
reg [2:0] total_error_count_w, total_error_count_r;

assign abs_llr[0] = (idata[63]) ? (~idata[63:56] + 1'b1) : idata[63:56]; // x^n
assign abs_llr[1] = (idata[55]) ? (~idata[55:48] + 1'b1) : idata[55:48];
assign abs_llr[2] = (idata[47]) ? (~idata[47:40] + 1'b1) : idata[47:40];
assign abs_llr[3] = (idata[39]) ? (~idata[39:32] + 1'b1) : idata[39:32];
assign abs_llr[4] = (idata[31]) ? (~idata[31:24] + 1'b1) : idata[31:24];
assign abs_llr[5] = (idata[23]) ? (~idata[23:16] + 1'b1) : idata[23:16];
assign abs_llr[6] = (idata[15]) ? (~idata[15:8] + 1'b1) : idata[15:8];
assign abs_llr[7] = (idata[7]) ? (~idata[7:0] + 1'b1) : idata[7:0]; //x^(n-7)

assign syndrome_out[3][0] = syndrome_out[0][0] ^ syndrome_out[1][0] ^ syndrome_out[2][0];
assign syndrome_out[3][1] = syndrome_out[0][1] ^ syndrome_out[1][1] ^ syndrome_out[2][1];
assign syndrome_out[3][2] = syndrome_out[0][2] ^ syndrome_out[1][2] ^ syndrome_out[2][2];
assign syndrome_out[3][3] = syndrome_out[0][3] ^ syndrome_out[1][3] ^ syndrome_out[2][3];
assign syndrome_out[3][4] = syndrome_out[0][4] ^ syndrome_out[1][4] ^ syndrome_out[2][4];
assign syndrome_out[3][5] = syndrome_out[0][5] ^ syndrome_out[1][5] ^ syndrome_out[2][5];
assign syndrome_out[3][6] = syndrome_out[0][6] ^ syndrome_out[1][6] ^ syndrome_out[2][6];
assign syndrome_out[3][7] = syndrome_out[0][7] ^ syndrome_out[1][7] ^ syndrome_out[2][7];

integer i, j, k;
// state machine 
always @(*) begin
	counter_w = counter_r;
	case(state_r)
		IDLE: begin
			if(set) begin
				state_w = SYNDROME_GEN;
				counter_w = 10'd0;
			end
			else begin
				state_w = IDLE;
				counter_w = 10'd0;
			end
		end
		SYNDROME_GEN: begin
			if(ready_r) begin
				counter_w = counter_r + 1;
			end
			else if(counter_r >= input_done_r -3) begin
				counter_w = counter_r + 1;
			end else begin
				counter_w = counter_r;
			end
			if(counter_r == input_done_r) begin
				state_w = ERROR_POLY_DISCREPENCY;
				counter_w = 0;
			end
			else begin
				state_w = SYNDROME_GEN;
			end
		end
		ERROR_POLY_DISCREPENCY: begin
			counter_w = counter_r + 1;
			if((counter_r == 0)) begin
				if(code_r == 2'b11) begin
					// early stop for BCH(1023,983)
					if((current_syndrome_r[0][0] == 10'd0) && (syndrome_r[0][0] == 10'd0) && (syndrome_r[0][1] == 10'd0) && (syndrome_r[0][2] == 10'd0) && (syndrome_r[0][3] == 10'd0) &&
						(syndrome_r[0][4] == 10'd0) && (syndrome_r[0][5] == 10'd0) && (syndrome_r[0][6] == 10'd0)) begin
						state_w = OUTPUT;
						counter_w = 0;
					end else begin
						state_w = ERROR_POLY_UPDATE;
					end
				end else begin
					// early stop for BCH(255,239) and BCH(63,51)
					if((current_syndrome_r[0][0] == 10'd0) && (syndrome_r[0][0] == 10'd0) && (syndrome_r[0][1] == 10'd0) && (syndrome_r[0][2] == 10'd0)) begin
						state_w = OUTPUT;
						counter_w = 0;
					end else begin
						state_w = ERROR_POLY_UPDATE;
					end
				end
			end else if(code_r == 2'b11) begin
				if (counter_r == 8) begin // max 8 iterations for BCH(1023,983)
					state_w = ERROR_LOCATE;
					counter_w = 0;
				end else begin
					state_w = ERROR_POLY_UPDATE;
				end
			end else begin
				if (counter_r == 4) begin // max 4 iterations for BCH(255,239) and BCH(63,51)
					state_w = ERROR_LOCATE;
					counter_w = 0;
				end else begin
					state_w = ERROR_POLY_UPDATE;
				end
			end
		end
		ERROR_POLY_UPDATE: begin
			state_w = ERROR_POLY_DISCREPENCY;
		end
		ERROR_LOCATE: begin
			counter_w = counter_r + 1;
			if(counter_r == input_done_r-2) begin
				if(mode_r) begin
                    state_w = CORRELATION;
					counter_w = 0;
				end else begin
					state_w = OUTPUT;
					counter_w = 0;
				end
			end else if(!mode_r) begin
				case(code_r)
					2'b11: begin
						if(error_count_r[0] == 3'd4) begin
							state_w = OUTPUT;
							counter_w = 0;
						end else begin
							state_w = ERROR_LOCATE;
						end
					end
					default:begin
						if(error_count_r[0] == 3'd2) begin
							state_w = OUTPUT;
							counter_w = 0;
						end else begin
							state_w = ERROR_LOCATE;
						end
					end
				endcase
			end else begin
                state_w = ERROR_LOCATE;
            end
		end
		CORRELATION: begin
			//SOFT DECODING VIA LLR CORRELATION
			counter_w = counter_r + 1;
			if(counter_r > 1) begin
				state_w = OUTPUT;
                counter_w = 0;
			end else begin
				state_w = CORRELATION;
			end
		end
		OUTPUT: begin
            counter_w = counter_r + 1;
			if(counter_r == total_error_count_r-1 || total_error_count_r == 0) begin
                state_w = IDLE;
            end else begin
                state_w = OUTPUT;
            end
		end
		default: begin
			state_w = SYNDROME_GEN;
		end
	endcase
end

always @(*) begin
	// default assignments
	code_w = code_r;
	mode_w = mode_r;
    output_mode_w = output_mode_r;
	input_done_w = input_done_r;
	ready_w = ready_r;
	finish_w = 1'b0;
	syndrome_data_in_w = 8'd0;
    syndrome_data_in_piped_w = syndrome_data_in_r;
	min_idx0_w = min_idx0_r;
	min_idx1_w = min_idx1_r;
	min_val0_w = min_val0_r;
	min_val1_w = min_val1_r;
	init = 0;
	corr_index_w = corr_index_r;
	for(i = 0; i < 1024; i = i + 1) begin
		llr_w[i] = llr_r[i];
	end
	for(j = 0; j < 4; j = j + 1) begin
		for(i = 0; i < 7; i = i + 1) begin
			syndrome_w[j][i] = syndrome_r[j][i+1];
		end
		syndrome_w[j][7] = 10'd0;
	end
	for(i = 0; i < 3; i = i + 1) begin
		syndrome_data_flipped_w[i] = syndrome_data_in_piped_r;
	end

	for(k = 0; k < 4; k = k+1) begin
		for(i = 0; i < 5; i = i + 1) begin
			sigma_w[k][i] = sigma_r[k][i];
			current_syndrome_w[k][i] = current_syndrome_r[k][i];
			Bx_w[k][i] = Bx_r[k][i];
		end

		delta_w[k] = delta_r[k];
		gamma_w[k] = gamma_r[k];
		L_w[k] = L_r[k];
		m_w[k] = m_r[k];
		for(i = 0; i < 10; i = i + 1) begin
			gf_mult_in1[k][i] = 10'd0;
			gf_mult_in2[k][i] = 10'd0;
		end

		for(i = 0; i < 4; i = i + 1) begin
			error_position_w[k][i] = error_position_r[k][i];
		end
		error_count[k][6] = error_count_r[k];
		error_count[k][5] = error_count[k][6];
		error_count[k][4] = error_count[k][5];
		error_count[k][3] = error_count[k][4];
		error_count[k][2] = error_count[k][3];
		error_count[k][1] = error_count[k][2];
		error_count[k][0] = error_count[k][1];
		error_count_w[k] = error_count[k][0];

		for(i = 0; i < 4; i = i + 1) begin
			partial_sum_w[k][i] = partial_sum_r[k][i];
		end
        for(i = 0; i < 8; i = i + 1) begin
            partial_sum[k][i] = 0;
        end

		corr_w[k] = corr_r[k];
	end


	odata_w = odata_r;
    total_error_count_w = total_error_count_r;
	error_output_w = error_output_r;

	case(state_r)
		IDLE: begin
			if(set) begin
				code_w = code;
				mode_w = mode;
				ready_w = 1'b1;
                output_mode_w = mode;
				case(code)
					2'b01: input_done_w = 8'd10; // BCH(63,51) -> 64 bits / 8 bits per symbol = 8 symbols
					2'b10: input_done_w = 8'd34; // BCH(255,239) -> 256 bits / 8 bits per symbol = 32 symbols
					default: input_done_w = 8'd130; // BCH(1023,983) -> 1024 bits / 8 bits per symbol = 128 symbols
				endcase
			end else begin
				code_w = code_r;
				mode_w = mode_r;
				input_done_w = input_done_r;
			end
            for(k = 0; k < 4; k =k +1) begin
                for(i = 0; i < 4; i = i + 1) begin
                    error_position_w[k][i] = 10'd0;
                    error_count_w[k] = 3'd0;
                end
            end
			total_error_count_w = 3'd0;
		end
		SYNDROME_GEN: begin
			total_error_count_w = 3'd0;
			
			// shift register
			if(ready_r) begin
				for (i = 0; i < (1024-8); i = i + 1) begin
					llr_w[i+8] = llr_r[i];
				end
				llr_w[7] = abs_llr[0]; // x^n
				llr_w[6] = abs_llr[1];
				llr_w[5] = abs_llr[2];
				llr_w[4] = abs_llr[3];
				llr_w[3] = abs_llr[4];
				llr_w[2] = abs_llr[5];
				llr_w[1] = abs_llr[6];
				llr_w[0] = abs_llr[7]; // x^n-7
			end
			// LLR is in decending power 
			// LLR0 = idata[63:56] , LLR1 = idata[55:48] , ... , LLR7 = idata[7:0]
			// signed number representation -> + for 0 , - for 1 -> same as sign bit
			syndrome_data_in_w = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
            syndrome_data_in_piped_w = syndrome_data_in_r;
			for(i = 0; i < 8; i = i + 1) begin
				for(j = 0; j < 4; j = j + 1) begin
					syndrome_w[j][i] = syndrome_out[j][i];
				end
			end

			// find smallest 2 LLR magnitudes and their indices for soft decoding
			if(counter_r > 1 && counter_r < input_done_r) begin
				min_idx0_w = min_idx0;
				min_idx1_w = min_idx1;
				min_val0_w = min_val0;
				min_val1_w = min_val1;
			end else begin
				min_idx0_w = min_idx0_r;
				min_idx1_w = min_idx1_r;
				min_val0_w = min_val0_r;
				min_val1_w = min_val1_r;
			end

			syndrome_data_flipped_w[0] = syndrome_data_in_piped_r;

			// did the smallest index change
			if(counter_r >0) begin
				if(min_idx0 != min_idx0_r) begin
					syndrome_data_flipped_w[1] = syndrome_data_in_piped_r ^ (8'b1 << (3'd7- min_idx0[2:0]));
					for(i = 0; i < 8; i = i +1) begin
						syndrome_w[1][i] = syndrome_out[0][i];
					end
				end else begin
					syndrome_data_flipped_w[1] = syndrome_data_in_piped_r;
				end

				// did the second smallest index change
				if(min_idx1 != min_idx1_r) begin
					// is the second smallest the previous smallest
					if(min_idx1 == min_idx0_r) begin
						syndrome_data_flipped_w[2] = syndrome_data_in_piped_r;
						for(i = 0; i < 8; i = i +1) begin
							syndrome_w[2][i] = syndrome_out[1][i];
						end
					end else begin
						syndrome_data_flipped_w[2] = syndrome_data_in_piped_r ^ (8'b1 << (3'd7- min_idx1[2:0]));
						for(i = 0; i < 8; i = i +1) begin
							syndrome_w[2][i] = syndrome_out[0][i];
						end
					end
				end else begin
					syndrome_data_flipped_w[2] = syndrome_data_in_piped_r;
				end
			end

			if(counter_r >= input_done_r -3) begin
				ready_w = 0;
			end

			//initialize error locator polynomial: sigma(x) = 1
			for(k = 0; k < 4; k = k+1) begin
				sigma_w[k][0] = 10'd1;
				for(i = 1; i < 5; i = i + 1) begin
					sigma_w[k][i] = 10'd0;
				end
				// initialize b(x) = 1
				Bx_w[k][0] = 10'd1;
				for(i = 1; i < 5; i = i + 1) begin
					Bx_w[k][i] = 10'd0;
				end
				// initialize other variables
				L_w[k] = 4'd0;
				m_w[k] = 4'd0;
				delta_w[k] = 10'd0;
				gamma_w[k] = 10'd1;
				
				// initialize current syndrome and shift syndrome registers
				if(counter_r == input_done_r ) begin
					current_syndrome_w[k][0] = syndrome_out[k][0];
					for(i = 0; i < 4; i = i + 1) begin
						current_syndrome_w[k][i+1] = current_syndrome_r[k][i];
					end
					for(i = 0; i < 7; i = i + 1) begin
						syndrome_w[k][i] = syndrome_out[k][i+1];
					end
					syndrome_w[k][7] = 10'd0;
				end
			end
		end
		ERROR_POLY_DISCREPENCY: begin
			total_error_count_w = 3'd0;
			// Inversionless Berlekamp-Massey algorithm to calculate error locator polynomial
			//calculate discrepancy
			if((counter_r == 0)) begin
				if(code_r == 2'b11) begin
					// early stop for BCH(1023,983)
					if((current_syndrome_r[0][0] == 10'd0) && (syndrome_r[0][0] == 10'd0) && (syndrome_r[0][1] == 10'd0) && (syndrome_r[0][2] == 10'd0) && (syndrome_r[0][3] == 10'd0) &&
						(syndrome_r[0][4] == 10'd0) && (syndrome_r[0][5] == 10'd0) && (syndrome_r[0][6] == 10'd0)) begin
						output_mode_w = 1'b0;
					end
				end else begin
					// early stop for BCH(255,239) and BCH(63,51)
					if((current_syndrome_r[0][0] == 10'd0) && (syndrome_r[0][0] == 10'd0) && (syndrome_r[0][1] == 10'd0) && (syndrome_r[0][2] == 10'd0)) begin
						output_mode_w = 1'b0;
					end
				end
			end
			for(k = 0; k < 4; k = k+1) begin
				for(i = 0; i < 5; i = i + 1) begin 
					gf_mult_in1[k][i] = sigma_r[k][i];
					gf_mult_in2[k][i] = current_syndrome_r[k][i];
				end

				delta_w[k] = gf_mult_out[k][0] ^ gf_mult_out[k][1] ^ gf_mult_out[k][2] ^ gf_mult_out[k][3] ^ gf_mult_out[k][4];

				// update current syndrome and shift syndrome registers
				current_syndrome_w[k][0] = syndrome_r[k][0];
				for(i = 0; i < 4; i = i + 1) begin
					current_syndrome_w[k][i+1] = current_syndrome_r[k][i];
				end
				for(i = 0; i < 7; i = i + 1) begin
					syndrome_w[k][i] = syndrome_r[k][i+1];
				end
				syndrome_w[k][7] = 10'd0;
				// update iteration counter
				m_w[k] = m_r[k] + 1;
			end
		end
		ERROR_POLY_UPDATE: begin
			// update scaling factor, degree of error locator polynomial, previous error locator polynomial
			// update error locator polynomial
			// sigma(x) = gamma*sigma(x) - delta * x * B(x)
			for(k = 0; k < 4; k = k+1) begin
				for(i = 0; i < 8; i = i + 1) begin
					syndrome_w[k][i] = syndrome_r[k][i];
				end
				for(i = 0; i < 5; i = i + 1) begin
					gf_mult_in1[k][i] = gamma_r[k];
					gf_mult_in2[k][i] = sigma_r[k][i];
				end
				for(i = 0; i < 5; i = i + 1) begin
					gf_mult_in1[k][i+5] = delta_r[k];
					gf_mult_in2[k][i+5] = Bx_r[k][i];
				end
				// calculate sigma(x)
				for(i = 1; i < 5; i = i + 1) begin
					sigma_w[k][i] = gf_mult_out[k][i] ^ gf_mult_out[k][(i-1)+5]; // x*B(x) -> shift by 1 -> (i-1)
				end
				sigma_w[k][0] = gf_mult_out[k][0]; // no x*B(x) term for x^0

				// update other variables
				if(delta_r[k] != 10'd0 && (L_r[k] << 1) <= m_r[k]) begin
					// if delta != 0 and 2*L <= m; 
					// B = sigma
					// L = m + 1 - L
					// gamma = delta
					for(i = 0; i < 5; i = i + 1) begin
						Bx_w[k][i] = sigma_r[k][i];
					end
					L_w[k] = m_r[k] - L_r[k];
					gamma_w[k] = delta_r[k];
				end else begin
					// else;
					// B = x*B
					// L, gamma unchanged
					for(i = 4; i > 0; i = i - 1) begin
						Bx_w[k][i] = Bx_r[k][i-1];
					end
					Bx_w[k][0] = 10'd0;
					L_w[k] = L_r[k];
					gamma_w[k] = gamma_r[k];
				end
			end
		end
		ERROR_LOCATE: begin
            output_mode_w = mode_r;
			// Chien Search to locate error positions
			// 2 parallel computation
			// 8 parallel iterations -> initialize inputs for first computation
			// first computation:
			// sigma(alpha^1) => position n - 1
			// sigma(alpha^2) => position n - 2
			// sigma(alpha^3) => position n - 3 //shift 3, 6, 9, 12
			// sigma(alpha^4) => position n - 4 //shift 4, 8, 12, 16
			// sigma(alpha^5) => position n - 5
			// sigma(alpha^6) => position n - 6
			// sigma(alpha^7) => position n - 7
			// sigma(alpha^8) => position n - 8

			// sigma(alpha^1 * alpha^(8* counter)) => position n - 1 - 8*counter     sigma == sigma(alpha^1)
			// sigma(alpha^2 * alpha^(8* counter)) => position n - 2 - 8*counter

			// a^0 = a^-n = pos(0), a^1 = a^( - (n-1)) = pos(n-1)
			// only valid if counter > 0;
			// counter = 1:
			// error[i] -> n - 1 - i
			// counter = k:
			// error[i] -> n - (1 + 8*(k-1) + i) = n - 1 - 8*(k-1) - i = n - 8*k + 7 - i

			
			
			if(counter_r > 0) begin
				for (i = 0; i < (1024-8); i = i + 1) begin
					llr_w[i+8] = llr_r[i];
				end
				for(k = 0; k < 4; k = k+1) begin
                    if(error_found_r[k][7]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count_r[k]] = 11'd70 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][7] = llr_r[62]; // 62, 54, 46, 38, 30, 22, 14, 6
							end
							2'b10: begin 
								error_position_w[k][error_count_r[k]] = 11'd262 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][7] = llr_r[254];
							end
							default: begin
								error_position_w[k][error_count_r[k]] = 11'd1030 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][7] = llr_r[1022];
							end
						endcase
						error_count[k][6] = error_count_r[k] + 3'd1;
					end else begin
						error_count[k][6] = error_count_r[k];
						partial_sum[k][7] = 7'd0;
					end
                    if(error_found_r[k][6]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][6]] = 11'd69 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][6] = llr_r[61];
							end
							2'b10: begin
								error_position_w[k][error_count[k][6]] = 11'd261 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][6] = llr_r[253];
							end
							default: begin
								error_position_w[k][error_count[k][6]] = 11'd1029 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][6] = llr_r[1021];
							end
						endcase
						error_count[k][5] = error_count[k][6] + 3'd1;
					end else begin
						error_count[k][5] = error_count[k][6];
						partial_sum[k][6] = 7'd0;
					end
                    if(error_found_r[k][5]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][5]] = 11'd68 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][5] = llr_r[60];
							end
							2'b10: begin
								error_position_w[k][error_count[k][5]] = 11'd260 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][5] = llr_r[252];
							end
							default: begin
								error_position_w[k][error_count[k][5]] = 11'd1028 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][5] = llr_r[1020];
							end
						endcase
						error_count[k][4] = error_count[k][5] + 3'd1;
					end else begin
						error_count[k][4] = error_count[k][5];
						partial_sum[k][5] = 7'd0;
					end
                    if(error_found_r[k][4]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][4]] = 11'd67 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][4] = llr_r[59];
							end
							2'b10: begin
								error_position_w[k][error_count[k][4]] = 11'd259 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][4] = llr_r[251];
							end
							default: begin
								error_position_w[k][error_count[k][4]] = 11'd1027 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][4] = llr_r[1019];
							end
						endcase
						error_count[k][3] = error_count[k][4] + 3'd1;
					end else begin
						error_count[k][3] = error_count[k][4];
						partial_sum[k][4] = 7'd0;
					end
                    if(error_found_r[k][3]) begin
						case(code_r)
							2'b01: begin 
								error_position_w[k][error_count[k][3]] = 11'd66 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][3] = llr_r[58];
							end
							2'b10: begin
								error_position_w[k][error_count[k][3]] = 11'd258 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][3] = llr_r[250];
							end
							default: begin
								error_position_w[k][error_count[k][3]] = 11'd1026 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][3] = llr_r[1018];
							end
						endcase
						error_count[k][2] = error_count[k][3] + 3'd1;
					end else begin
						error_count[k][2] = error_count[k][3];
						partial_sum[k][3] = 7'd0;
					end
                    if(error_found_r[k][2]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][2]] = 11'd65 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][2] = llr_r[57];
							end
							2'b10: begin
								error_position_w[k][error_count[k][2]] = 11'd257 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][2] = llr_r[249];
							end
							default: begin
								error_position_w[k][error_count[k][2]] = 11'd1025 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][2] = llr_r[1017];
							end
						endcase
						error_count[k][1] = error_count[k][2] + 3'd1;
					end else begin
						error_count[k][1] = error_count[k][2];
						partial_sum[k][2] = 7'd0;
					end
                    if(error_found_r[k][1]) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][1]] = 11'd64 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][1] = llr_r[56];
							end
							2'b10: begin
								error_position_w[k][error_count[k][1]] = 11'd256 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][1] = llr_r[248];
							end
							default: begin
								error_position_w[k][error_count[k][1]] = 11'd1024 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][1] = llr_r[1016];
							end
						endcase
						error_count[k][0] = error_count[k][1] + 3'd1;
					end else begin
						error_count[k][0] = error_count[k][1];
						partial_sum[k][1] = 7'd0;
					end
					// ignore error_found[k][0] at last cycle becuase it is -1 position
					if(error_found_r[k][0] && (counter_r != input_done_r-2) ) begin
						case(code_r)
							2'b01: begin
								error_position_w[k][error_count[k][0]] = 11'd63 - (counter_r << 3); // BCH(63,51)
								partial_sum[k][0] = llr_r[55];
							end
							2'b10: begin
								error_position_w[k][error_count[k][0]] = 11'd255 - (counter_r << 3); // BCH(255,239)
								partial_sum[k][0] = llr_r[247];
							end
							default: begin
								error_position_w[k][error_count[k][0]] = 11'd1023 - (counter_r << 3); // BCH(1023,983)
								partial_sum[k][0] = llr_r[1015];
							end
						endcase
						error_count_w[k] = error_count[k][0] + 3'd1;
					end else begin
						error_count_w[k] = error_count[k][0];
						partial_sum[k][0] = 7'd0;
					end
				end
			end else begin
				partial_sum[1][0] = min_val0_r;
				partial_sum[2][0] = min_val1_r;
				partial_sum[3][0] = min_val0_r;
				partial_sum[3][1] = min_val1_r;
				case(code_r)
					2'b01: begin // BCH(63,51)
						min_idx0_w = 10'd63-min_idx0_r;
						min_idx1_w = 10'd63-min_idx1_r;
					end
					2'b10: begin // BCH(255,239)
						min_idx0_w = 10'd255-min_idx0_r;
						min_idx1_w = 10'd255-min_idx1_r;
					end
					default: begin // BCH(1023,983)
						min_idx0_w = 10'd1023-min_idx0_r;
						min_idx1_w = 10'd1023-min_idx1_r;
					end
				endcase
			end

            for(k = 0; k < 4; k = k+1) begin
				corr_w[k] = corr_r[k] + partial_sum_r[k][0]+ partial_sum_r[k][1] + partial_sum_r[k][2] + partial_sum_r[k][3];
                partial_sum_w[k][0] = partial_sum[k][0] + partial_sum[k][1] ;
                partial_sum_w[k][1] = partial_sum[k][2] + partial_sum[k][3] ;
                partial_sum_w[k][2] = partial_sum[k][4] + partial_sum[k][5] ;
                partial_sum_w[k][3] = partial_sum[k][6] + partial_sum[k][7] ;
            end

			for(k = 0; k < 4; k = k+1) begin
				//update sigma_r for next iteration
				for(i = 1; i < 5; i = i + 1) begin
					sigma_w[k][i] = chien_shifted_value[k][7][i];
				end
				sigma_w[k][0] = sigma_r[k][0];
			end
            total_error_count_w = error_count_w[0];
		end
		CORRELATION: begin
			// choose the error based on LLR correlation
			corr_index_w = corr_index;
			if(counter_r > 1) begin
                case(corr_index)
                    2'b00: begin
                        total_error_count_w = error_count_r[corr_index_r];
                    end
                    2'b01: begin
                        total_error_count_w = error_count_r[corr_index_r]+1;
                    end
                    2'b10: begin
                        total_error_count_w = error_count_r[corr_index_r]+1;
                    end
                    default: begin
                        total_error_count_w = error_count_r[corr_index_r]+2;
                    end
                endcase
				error_output_w = error_position_r[corr_index_r][fetch_error_pos];
				error_count_w[corr_index_r] = error_count_r[corr_index_r] - 3'd1;
				if(error_count_r[0] == 3'd0) begin
					output_mode_w = 1'b0; // switch to hard decoding mode
                    error_count_w[0] = 0;
                    // error_output_w = 10'd1023; // no more errors to output
				end
			end else if (counter_r == 0) begin
				for(k = 0; k < 4; k = k+1) begin
					corr_w[k] = corr_r[k] + partial_sum_r[k][0]+ partial_sum_r[k][1] + partial_sum_r[k][2] + partial_sum_r[k][3];
                    partial_sum_w[k][0] = partial_sum[k][0] + partial_sum[k][1] ;
                    partial_sum_w[k][1] = partial_sum[k][2] + partial_sum[k][3] ;
                    partial_sum_w[k][2] = partial_sum[k][4] + partial_sum[k][5] ;
                    partial_sum_w[k][3] = partial_sum[k][6] + partial_sum[k][7] ;
                end
				if(error_count_r[0] != L_r[0] ) begin
					corr_w[0] = 10'd1023;
					error_count_w[0] = 1;
				end
				if(error_count_r[1] != L_r[1] ) begin
					corr_w[1] = 10'd1023;
					error_count_w[1] = 1;
				end
				if(error_count_r[2] != L_r[2] ) begin
					corr_w[2] = 10'd1023;
					error_count_w[2] = 1;
				end
				if(error_count_r[3] != L_r[3]) begin
					corr_w[3] = 10'd1023;
					error_count_w[3] = 1;
				end
			end

		end
		OUTPUT: begin
			// output corrected data
			error_output_w = error_output_r;
			finish_w = 1'b1;
			if(output_mode_r) begin
				case(corr_index_r)
					2'b00: begin
						odata_w = error_output_w;
                        error_count_w[0] = error_count_r[0] - 3'd1;
						error_output_w = error_position_r[corr_index_r][fetch_error_pos];
					end
					2'b01: begin
						if(error_output_r < min_idx0_r && (error_count_r[1] != 7)) begin
							odata_w = error_output_r; // if same as first index, choose first index
							error_count_w[1] = error_count_r[1] - 3'd1;
							error_output_w = error_position_r[corr_index_r][fetch_error_pos];
						end else begin
							odata_w = min_idx0_r;
							min_idx0_w = 10'd1023;
						end
					end
					2'b10: begin
						if(error_output_r < min_idx1_r && (error_count_r[2] != 7)) begin
							odata_w = error_output_r; // if same as second index, choose first index
							error_count_w[2] = error_count_r[2] - 3'd1;
							error_output_w = error_position_r[corr_index_r][fetch_error_pos];
						end else begin
							odata_w = min_idx1_r;
							min_idx1_w = 10'd1023;
						end
					end
					default: begin
						if(error_output_r < min_idx0_r && error_output_r < min_idx1_r && (error_count_r[3] != 7)) begin
							odata_w = error_output_r; // if same as first or second index, choose first index
							error_count_w[3] = error_count_r[3] - 3'd1;
							error_output_w = error_position_r[corr_index_r][fetch_error_pos];
						end else if(min_idx0_r < min_idx1_r) begin
							odata_w = min_idx0_r;
							min_idx0_w = 10'd1023;
						end else begin
							odata_w = min_idx1_r;
							min_idx1_w = 10'd1023;
						end
					end
				endcase
			end else begin
				error_count_w[0] = error_count_r[0] - 3'd1;
				if(error_count_r[0] == 3'd0) begin
					odata_w = 10'd1023; // no more errors to output
					min_idx0_w = 10'd1023;
					min_idx1_w = 10'd1023;
					min_val0_w = 7'd127;
					min_val1_w = 7'd127;
				end else begin
					odata_w = error_position_r[0][fetch_error_pos]; // hard decoding -> just output error positions
				end
				
			end

			// reset registers for next codeword
			// for(k = 0; k < 4; k = k+1) begin
			// 	for(i = 0; i < 8; i = i + 1) begin
			// 		syndrome_w[k][i] = 10'd0;
			// 	end
			// end
			// for(i = 0; i < 3; i = i + 1) begin
			// 	syndrome_data_flipped_w[i] = 8'd0;
			// end
            syndrome_data_in_piped_w = 8'd0;
			if(counter_r == total_error_count_r-1) begin          
				min_idx0_w = 10'd1023;
				min_idx1_w = 10'd1023;
				min_val0_w = 7'd127;
				min_val1_w = 7'd127;
			end
			for(k = 0; k < 4; k = k+1) begin
            	corr_w[k] = 10'd0;
			end
		end
	endcase
end

SYNDROME_GEN SYNDROME_GEN_00_INST(
	.i_syndrome_1(syndrome_r[0][0]),
	.i_syndrome_2(syndrome_r[0][1]),
	.i_syndrome_3(syndrome_r[0][2]),
	.i_syndrome_4(syndrome_r[0][3]),
	.i_syndrome_5(syndrome_r[0][4]),
	.i_syndrome_6(syndrome_r[0][5]),
	.i_syndrome_7(syndrome_r[0][6]),
	.i_syndrome_8(syndrome_r[0][7]),
	.syndrome_data_in(syndrome_data_flipped_r[0]),
	.code(code_r),
	.o_syndrome_1(syndrome_out[0][0]),
	.o_syndrome_2(syndrome_out[0][1]),
	.o_syndrome_3(syndrome_out[0][2]),
	.o_syndrome_4(syndrome_out[0][3]),
	.o_syndrome_5(syndrome_out[0][4]),
	.o_syndrome_6(syndrome_out[0][5]),
	.o_syndrome_7(syndrome_out[0][6]),
	.o_syndrome_8(syndrome_out[0][7])
);

SYNDROME_GEN SYNDROME_GEN_01_INST(
	.i_syndrome_1(syndrome_r[1][0]),
	.i_syndrome_2(syndrome_r[1][1]),
	.i_syndrome_3(syndrome_r[1][2]),
	.i_syndrome_4(syndrome_r[1][3]),
	.i_syndrome_5(syndrome_r[1][4]),
	.i_syndrome_6(syndrome_r[1][5]),
	.i_syndrome_7(syndrome_r[1][6]),
	.i_syndrome_8(syndrome_r[1][7]),
	.syndrome_data_in(syndrome_data_flipped_r[1]),
	.code(code_r),
	.o_syndrome_1(syndrome_out[1][0]),
	.o_syndrome_2(syndrome_out[1][1]),
	.o_syndrome_3(syndrome_out[1][2]),
	.o_syndrome_4(syndrome_out[1][3]),
	.o_syndrome_5(syndrome_out[1][4]),
	.o_syndrome_6(syndrome_out[1][5]),
	.o_syndrome_7(syndrome_out[1][6]),
	.o_syndrome_8(syndrome_out[1][7])
);

SYNDROME_GEN SYNDROME_GEN_10_INST(
	.i_syndrome_1(syndrome_r[2][0]),
	.i_syndrome_2(syndrome_r[2][1]),
	.i_syndrome_3(syndrome_r[2][2]),
	.i_syndrome_4(syndrome_r[2][3]),
	.i_syndrome_5(syndrome_r[2][4]),
	.i_syndrome_6(syndrome_r[2][5]),
	.i_syndrome_7(syndrome_r[2][6]),
	.i_syndrome_8(syndrome_r[2][7]),
	.syndrome_data_in(syndrome_data_flipped_r[2]),
	.code(code_r),
	.o_syndrome_1(syndrome_out[2][0]),
	.o_syndrome_2(syndrome_out[2][1]),
	.o_syndrome_3(syndrome_out[2][2]),
	.o_syndrome_4(syndrome_out[2][3]),
	.o_syndrome_5(syndrome_out[2][4]),
	.o_syndrome_6(syndrome_out[2][5]),
	.o_syndrome_7(syndrome_out[2][6]),
	.o_syndrome_8(syndrome_out[2][7])
);

BERLEKAMP_MULT BERLEKAMP_MULT00(
	.i_gf_mult_x1_0(gf_mult_in1[0][0]),
	.i_gf_mult_x1_1(gf_mult_in1[0][1]),
	.i_gf_mult_x1_2(gf_mult_in1[0][2]),
	.i_gf_mult_x1_3(gf_mult_in1[0][3]),
	.i_gf_mult_x1_4(gf_mult_in1[0][4]),
	.i_gf_mult_x1_5(gf_mult_in1[0][5]),
	.i_gf_mult_x1_6(gf_mult_in1[0][6]),
	.i_gf_mult_x1_7(gf_mult_in1[0][7]),
	.i_gf_mult_x1_8(gf_mult_in1[0][8]),
	.i_gf_mult_x1_9(gf_mult_in1[0][9]),

	.i_gf_mult_x2_0(gf_mult_in2[0][0]),
	.i_gf_mult_x2_1(gf_mult_in2[0][1]),
	.i_gf_mult_x2_2(gf_mult_in2[0][2]),
	.i_gf_mult_x2_3(gf_mult_in2[0][3]),
	.i_gf_mult_x2_4(gf_mult_in2[0][4]),
	.i_gf_mult_x2_5(gf_mult_in2[0][5]),
	.i_gf_mult_x2_6(gf_mult_in2[0][6]),
	.i_gf_mult_x2_7(gf_mult_in2[0][7]),
	.i_gf_mult_x2_8(gf_mult_in2[0][8]),
	.i_gf_mult_x2_9(gf_mult_in2[0][9]),
	.code(code_r),

	.o_gf_mult_y0(gf_mult_out[0][0]),
	.o_gf_mult_y1(gf_mult_out[0][1]),
	.o_gf_mult_y2(gf_mult_out[0][2]),
	.o_gf_mult_y3(gf_mult_out[0][3]),
	.o_gf_mult_y4(gf_mult_out[0][4]),
	.o_gf_mult_y5(gf_mult_out[0][5]),
	.o_gf_mult_y6(gf_mult_out[0][6]),
	.o_gf_mult_y7(gf_mult_out[0][7]),
	.o_gf_mult_y8(gf_mult_out[0][8]),
	.o_gf_mult_y9(gf_mult_out[0][9])
);

BERLEKAMP_MULT BERLEKAMP_MULT01(
	.i_gf_mult_x1_0(gf_mult_in1[1][0]),
	.i_gf_mult_x1_1(gf_mult_in1[1][1]),
	.i_gf_mult_x1_2(gf_mult_in1[1][2]),
	.i_gf_mult_x1_3(gf_mult_in1[1][3]),
	.i_gf_mult_x1_4(gf_mult_in1[1][4]),
	.i_gf_mult_x1_5(gf_mult_in1[1][5]),
	.i_gf_mult_x1_6(gf_mult_in1[1][6]),
	.i_gf_mult_x1_7(gf_mult_in1[1][7]),
	.i_gf_mult_x1_8(gf_mult_in1[1][8]),
	.i_gf_mult_x1_9(gf_mult_in1[1][9]),

	.i_gf_mult_x2_0(gf_mult_in2[1][0]),
	.i_gf_mult_x2_1(gf_mult_in2[1][1]),
	.i_gf_mult_x2_2(gf_mult_in2[1][2]),
	.i_gf_mult_x2_3(gf_mult_in2[1][3]),
	.i_gf_mult_x2_4(gf_mult_in2[1][4]),
	.i_gf_mult_x2_5(gf_mult_in2[1][5]),
	.i_gf_mult_x2_6(gf_mult_in2[1][6]),
	.i_gf_mult_x2_7(gf_mult_in2[1][7]),
	.i_gf_mult_x2_8(gf_mult_in2[1][8]),
	.i_gf_mult_x2_9(gf_mult_in2[1][9]),
	.code(code_r),

	.o_gf_mult_y0(gf_mult_out[1][0]),
	.o_gf_mult_y1(gf_mult_out[1][1]),
	.o_gf_mult_y2(gf_mult_out[1][2]),
	.o_gf_mult_y3(gf_mult_out[1][3]),
	.o_gf_mult_y4(gf_mult_out[1][4]),
	.o_gf_mult_y5(gf_mult_out[1][5]),
	.o_gf_mult_y6(gf_mult_out[1][6]),
	.o_gf_mult_y7(gf_mult_out[1][7]),
	.o_gf_mult_y8(gf_mult_out[1][8]),
	.o_gf_mult_y9(gf_mult_out[1][9])
);

BERLEKAMP_MULT BERLEKAMP_MULT10(
	.i_gf_mult_x1_0(gf_mult_in1[2][0]),
	.i_gf_mult_x1_1(gf_mult_in1[2][1]),
	.i_gf_mult_x1_2(gf_mult_in1[2][2]),
	.i_gf_mult_x1_3(gf_mult_in1[2][3]),
	.i_gf_mult_x1_4(gf_mult_in1[2][4]),
	.i_gf_mult_x1_5(gf_mult_in1[2][5]),
	.i_gf_mult_x1_6(gf_mult_in1[2][6]),
	.i_gf_mult_x1_7(gf_mult_in1[2][7]),
	.i_gf_mult_x1_8(gf_mult_in1[2][8]),
	.i_gf_mult_x1_9(gf_mult_in1[2][9]),
	.i_gf_mult_x2_0(gf_mult_in2[2][0]),
	.i_gf_mult_x2_1(gf_mult_in2[2][1]),
	.i_gf_mult_x2_2(gf_mult_in2[2][2]),
	.i_gf_mult_x2_3(gf_mult_in2[2][3]),
	.i_gf_mult_x2_4(gf_mult_in2[2][4]),
	.i_gf_mult_x2_5(gf_mult_in2[2][5]),
	.i_gf_mult_x2_6(gf_mult_in2[2][6]),
	.i_gf_mult_x2_7(gf_mult_in2[2][7]),
	.i_gf_mult_x2_8(gf_mult_in2[2][8]),
	.i_gf_mult_x2_9(gf_mult_in2[2][9]),

	.code(code_r),

	.o_gf_mult_y0(gf_mult_out[2][0]),
	.o_gf_mult_y1(gf_mult_out[2][1]),
	.o_gf_mult_y2(gf_mult_out[2][2]),
	.o_gf_mult_y3(gf_mult_out[2][3]),
	.o_gf_mult_y4(gf_mult_out[2][4]),
	.o_gf_mult_y5(gf_mult_out[2][5]),
	.o_gf_mult_y6(gf_mult_out[2][6]),
	.o_gf_mult_y7(gf_mult_out[2][7]),
	.o_gf_mult_y8(gf_mult_out[2][8]),
	.o_gf_mult_y9(gf_mult_out[2][9])
);



BERLEKAMP_MULT BERLEKAMP_MULT11(
	.i_gf_mult_x1_0(gf_mult_in1[3][0]),
	.i_gf_mult_x1_1(gf_mult_in1[3][1]),
	.i_gf_mult_x1_2(gf_mult_in1[3][2]),
	.i_gf_mult_x1_3(gf_mult_in1[3][3]),
	.i_gf_mult_x1_4(gf_mult_in1[3][4]),
	.i_gf_mult_x1_5(gf_mult_in1[3][5]),
	.i_gf_mult_x1_6(gf_mult_in1[3][6]),
	.i_gf_mult_x1_7(gf_mult_in1[3][7]),
	.i_gf_mult_x1_8(gf_mult_in1[3][8]),
	.i_gf_mult_x1_9(gf_mult_in1[3][9]),
	.i_gf_mult_x2_0(gf_mult_in2[3][0]),
	.i_gf_mult_x2_1(gf_mult_in2[3][1]),
	.i_gf_mult_x2_2(gf_mult_in2[3][2]),
	.i_gf_mult_x2_3(gf_mult_in2[3][3]),
	.i_gf_mult_x2_4(gf_mult_in2[3][4]),
	.i_gf_mult_x2_5(gf_mult_in2[3][5]),
	.i_gf_mult_x2_6(gf_mult_in2[3][6]),
	.i_gf_mult_x2_7(gf_mult_in2[3][7]),
	.i_gf_mult_x2_8(gf_mult_in2[3][8]),
	.i_gf_mult_x2_9(gf_mult_in2[3][9]),

	.code(code_r),

	.o_gf_mult_y0(gf_mult_out[3][0]),
	.o_gf_mult_y1(gf_mult_out[3][1]),
	.o_gf_mult_y2(gf_mult_out[3][2]),
	.o_gf_mult_y3(gf_mult_out[3][3]),
	.o_gf_mult_y4(gf_mult_out[3][4]),
	.o_gf_mult_y5(gf_mult_out[3][5]),
	.o_gf_mult_y6(gf_mult_out[3][6]),
	.o_gf_mult_y7(gf_mult_out[3][7]),
	.o_gf_mult_y8(gf_mult_out[3][8]),
	.o_gf_mult_y9(gf_mult_out[3][9])
);


genvar g;
generate
	for(g = 0; g < 4; g = g + 1) begin : gf_mul_alpha_1_gen

	// first row, multiply by alpha^k, k = 1-8
		gf_mul_alpha_1 gf_mul_alpha_1_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][0][1])
		);
		gf_mul_alpha_2 gf_mul_alpha_2_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][1][1])
		);
		gf_mul_alpha_3 gf_mul_alpha_3_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][2][1])
		);
		gf_mul_alpha_4 gf_mul_alpha_4_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][3][1])
		);
		gf_mul_alpha_5 gf_mul_alpha_5_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][4][1])
		);
		gf_mul_alpha_6 gf_mul_alpha_6_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][5][1])
		);
		gf_mul_alpha_7 gf_mul_alpha_7_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][6][1])
		);
		gf_mul_alpha_8 gf_mul_alpha_8_INST (
			.current(sigma_r[g][1]),
			.m(code_r),
			.next(chien_shifted_value[g][7][1])
		);
		// second row, multiply by alpha^(2k), k = 1-8
		gf_mul_alpha_2 gf_mul_alpha_2_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][0][2])
		);
		gf_mul_alpha_4 gf_mul_alpha_4_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][1][2])
		);
		gf_mul_alpha_6 gf_mul_alpha_6_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][2][2])
		);
		gf_mul_alpha_8 gf_mul_alpha_8_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][3][2])
		);
		gf_mul_alpha_10 gf_mul_alpha_10_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][4][2])
		);
		gf_mul_alpha_12 gf_mul_alpha_12_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][5][2])
		);
		gf_mul_alpha_14 gf_mul_alpha_14_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][6][2])
		);
		gf_mul_alpha_16 gf_mul_alpha_16_INST2 (
			.current(sigma_r[g][2]),
			.m(code_r),
			.next(chien_shifted_value[g][7][2])
		);
		// third row, multiply by alpha^(3k), k = 1-8
		gf_mul_alpha_3 gf_mul_alpha_3_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][0][3])
		);
		gf_mul_alpha_6 gf_mul_alpha_6_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][1][3])
		);
		gf_mul_alpha_9 gf_mul_alpha_9_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][2][3])
		);
		gf_mul_alpha_12 gf_mul_alpha_12_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][3][3])
		);
		gf_mul_alpha_15 gf_mul_alpha_15_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][4][3])
		);
		gf_mul_alpha_18 gf_mul_alpha_18_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][5][3])
		);
		gf_mul_alpha_21 gf_mul_alpha_21_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][6][3])
		);
		gf_mul_alpha_24 gf_mul_alpha_24_INST3 (
			.current(sigma_r[g][3]),
			.m(code_r),
			.next(chien_shifted_value[g][7][3])
		);
		// fourth row, multiply by alpha^(4k), k = 1-8
		gf_mul_alpha_4 gf_mul_alpha_4_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][0][4])
		);
		gf_mul_alpha_8 gf_mul_alpha_8_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][1][4])
		);
		gf_mul_alpha_12 gf_mul_alpha_12_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][2][4])
		);
		gf_mul_alpha_16 gf_mul_alpha_16_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][3][4])
		);
		gf_mul_alpha_20 gf_mul_alpha_20_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][4][4])
		);
		gf_mul_alpha_24 gf_mul_alpha_24_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][5][4])
		);
		gf_mul_alpha_28 gf_mul_alpha_28_INST4 (	
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][6][4])
		);
		gf_mul_alpha_32 gf_mul_alpha_32_INST4 (
			.current(sigma_r[g][4]),
			.m(code_r),
			.next(chien_shifted_value[g][7][4])
		);
	end
endgenerate	






MINCORR MINCORR_INST (
	.in_corr0(corr_r[0]),
	.in_corr1(corr_r[1]),
	.in_corr2(corr_r[2]),
	.in_corr3(corr_r[3]),
	.min_corr(corr_index)
);


MIN2_LLR MIN2_LLR_INST (
	// input [6:0] abs_llr [0:7],
	// input [6:0] i_min_val0,
	// input [6:0] i_min_val1,
	// input [9:0] i_min_idx0,
	// input [9:0] i_min_idx1,

    // input [9:0] counter,

    // output reg [9:0] o_min_idx0,
    // output reg [9:0] o_min_idx1,
	// output reg [6:0] o_min_val0,
	// output reg [6:0] o_min_val1
    .clk(clk),
    .rstn(rstn),
	.clear(clear),
	.abs_llr0((counter_r == 1) ? 7'd127 : llr_r[7]),
	.abs_llr1(llr_r[6]),
	.abs_llr2(llr_r[5]),
	.abs_llr3(llr_r[4]),
	.abs_llr4(llr_r[3]),
	.abs_llr5(llr_r[2]),
	.abs_llr6(llr_r[1]),
	.abs_llr7(llr_r[0]),
	.counter(counter_r),

	.i_min_val0(min_val0_r),
	.i_min_val1(min_val1_r),
	.i_min_idx0(min_idx0_r),
	.i_min_idx1(min_idx1_r),
	.o_min_idx0(min_idx0),
	.o_min_idx1(min_idx1),
	.o_min_val0(min_val0),
	.o_min_val1(min_val1)

);

wire clk_en = mode_r | ~rstn; // clock enable for soft decoding and reset
wire gclk = clk & clk_en; // gated clock for soft decoding

always @(posedge gclk) begin
    if(!rstn) begin
        for(k = 1; k < 4; k = k + 1) begin
            for(i = 0; i < 8; i = i + 1) begin
                syndrome_r[k][i] <= 10'd0;
            end
            for(i = 0; i < 5; i = i + 1) begin
                sigma_r[k][i] <= 10'd0;
                current_syndrome_r[k][i] <= 10'd0;
                Bx_r[k][i] <= 10'd0;
            end
            delta_r[k] <= 10'd0;
            gamma_r[k] <= 10'd1;
            L_r[k] <= 4'd0;
            m_r[k] <= 4'd0;
            for(i = 0; i < 4; i = i + 1) begin
                error_position_r[k][i] <= 10'd0;
            end
            error_count_r[k] <= 3'd0;
            corr_r[k] <= 10'd0;
            for(i = 0; i < 4; i = i + 1) begin
                partial_sum_r[k][i] <= 7'd0;
            end
            for(i = 0; i < 8; i = i + 1) begin
                error_found_r[k][i] <= 1'b0;
            end
		end
    end else begin
        for(k = 1; k < 4; k = k + 1) begin
            for(i = 0; i < 8; i = i + 1) begin
                syndrome_r[k][i] <= syndrome_w[k][i];
            end
            for(i = 0; i < 5; i = i + 1) begin
                sigma_r[k][i] <= sigma_w[k][i];
                current_syndrome_r[k][i] <= current_syndrome_w[k][i];
                Bx_r[k][i] <= Bx_w[k][i];
            end
            delta_r[k] <= delta_w[k];
            gamma_r[k] <= gamma_w[k];
            L_r[k] <= L_w[k];
            m_r[k] <= m_w[k];
            for(i = 0; i < 4; i = i + 1) begin
                error_position_r[k][i] <= error_position_w[k][i];
            end
            error_count_r[k] <= error_count_w[k];
            corr_r[k] <= corr_w[k];
            for(i = 0; i < 4; i = i + 1) begin
                partial_sum_r[k][i] <= partial_sum_w[k][i];
            end
            for(i = 0; i < 8; i = i + 1) begin
                error_found_r[k][i] <= error_found_w[k][i];
            end
		end
    end
end

// sequential logic
always @(posedge clk) begin
	if(!rstn) begin
		// reset all registers
		state_r <= IDLE;
		code_r <= 2'b00;
		mode_r <= 1'b0;
        output_mode_r <= 1'b0;
		counter_r <= 10'd0;
		ready_r <= 1'b0;
		finish_r <= 1'b0;
		syndrome_data_in_r <= 8'd0;
        syndrome_data_in_piped_r <= 8'd0;
		input_done_r <= 10'd0;
		for(i = 0; i < 1024; i = i + 1) begin
			llr_r[i] <= 8'd0;
		end

		for(k = 0; k < 3; k = k + 1) begin
			syndrome_data_flipped_r[k] <= 8'd0;
		end


        for(i = 0; i < 8; i = i + 1) begin
            syndrome_r[0][i] <= 10'd0;
        end

        for(i = 0; i < 5; i = i + 1) begin
            sigma_r[0][i] <= 10'd0;
            current_syndrome_r[0][i] <= 10'd0;
            Bx_r[0][i] <= 10'd0;
        end
        delta_r[0] <= 10'd0;
        gamma_r[0] <= 10'd1;
        L_r[0] <= 4'd0;
        m_r[0] <= 4'd0;
        for(i = 0; i < 4; i = i + 1) begin
            error_position_r[0][i] <= 10'd0;
        end
        error_count_r[0] <= 3'd0;
        corr_r[0] <= 10'd0;
        for(i = 0; i < 4; i = i + 1) begin
            partial_sum_r[0][i] <= 7'd0;
        end
        for(i = 0; i < 8; i = i + 1) begin
            error_found_r[0][i] <= 1'b0;
        end

		error_output_r <= 10'd0;
		corr_index_r <= 2'b00;
		min_idx0_r <= 10'd1023;
		min_idx1_r <= 10'd1023;
		min_val0_r <= 7'd127;
		min_val1_r <= 7'd127;

        total_error_count_r <= 10'd0;

		odata_r <= 10'd0;
	end else begin
		// update all registers
		state_r <= state_w;
		code_r <= code_w;
		mode_r <= mode_w;
        output_mode_r <= output_mode_w;
		counter_r <= counter_w;
		ready_r <= ready_w;
		finish_r <= finish_w;
		syndrome_data_in_r <= syndrome_data_in_w;
        syndrome_data_in_piped_r <= syndrome_data_in_piped_w;  
		input_done_r <= input_done_w;
		for(i = 0; i < 1024; i = i + 1) begin
			llr_r[i] <= llr_w[i];
		end

		for(k = 0; k < 3; k = k + 1) begin
			syndrome_data_flipped_r[k] <= syndrome_data_flipped_w[k];
		end

		// for(k = 1; k < 4; k = k + 1) begin
        //     if(mode_r) begin
        //         for(i = 0; i < 8; i = i + 1) begin
        //             syndrome_r[k][i] <= syndrome_w[k][i];
        //         end
        //         for(i = 0; i < 5; i = i + 1) begin
        //             sigma_r[k][i] <= sigma_w[k][i];
        //             current_syndrome_r[k][i] <= current_syndrome_w[k][i];
        //             Bx_r[k][i] <= Bx_w[k][i];
        //         end
        //         delta_r[k] <= delta_w[k];
        //         gamma_r[k] <= gamma_w[k];
        //         L_r[k] <= L_w[k];
        //         m_r[k] <= m_w[k];
        //         for(i = 0; i < 4; i = i + 1) begin
        //             error_position_r[k][i] <= error_position_w[k][i];
        //         end
        //         error_count_r[k] <= error_count_w[k];
        //         corr_r[k] <= corr_w[k];
        //         for(i = 0; i < 4; i = i + 1) begin
        //             partial_sum_r[k][i] <= partial_sum_w[k][i];
        //         end
        //         for(i = 0; i < 8; i = i + 1) begin
        //             error_found_r[k][i] <= error_found_w[k][i];
        //         end
        //     end
		// end
        for(i = 0; i < 8; i = i + 1) begin
            syndrome_r[0][i] <= syndrome_w[0][i];
        end
        for(i = 0; i < 5; i = i + 1) begin
            sigma_r[0][i] <= sigma_w[0][i];
            current_syndrome_r[0][i] <= current_syndrome_w[0][i];
            Bx_r[0][i] <= Bx_w[0][i];
        end
        delta_r[0] <= delta_w[0];
        gamma_r[0] <= gamma_w[0];
        L_r[0] <= L_w[0];
        m_r[0] <= m_w[0];
        for(i = 0; i < 4; i = i + 1) begin
            error_position_r[0][i] <= error_position_w[0][i];
        end
        error_count_r[0] <= error_count_w[0];
        corr_r[0] <= corr_w[0];
        for(i = 0; i < 4; i = i + 1) begin
            partial_sum_r[0][i] <= partial_sum_w[0][i];
        end
        for(i = 0; i < 8; i = i + 1) begin
            error_found_r[0][i] <= error_found_w[0][i];
        end
		error_output_r <= error_output_w;
		corr_index_r <= corr_index_w;
		min_idx0_r <= min_idx0_w;
		min_idx1_r <= min_idx1_w;
		min_val0_r <= min_val0_w;
		min_val1_r <= min_val1_w;
        total_error_count_r <= total_error_count_w;

		odata_r <= odata_w;
	end
end
endmodule


// Multiply a GF(2^m) element by alpha^8
// Supports m = 6, 8, 10
module chien_mul_a8 (
    input  [9:0] v_in,   // s0..s9 (only lower m bits used)
    input  [1:0] m,      // 6, 8, or 10
    output [9:0] v_out
);
    // Alias bits for readability
    wire s0 = v_in[0];
    wire s1 = v_in[1];
    wire s2 = v_in[2];
    wire s3 = v_in[3];
    wire s4 = v_in[4];
    wire s5 = v_in[5];
    wire s6 = v_in[6];
    wire s7 = v_in[7];
    wire s8 = v_in[8];
    wire s9 = v_in[9];

    reg [9:0] next;

    always @* begin
        case (m)
            2'b01: begin
                // m = 6, p(x) = x^6 + x + 1, multiply by alpha^8
                next[0] = s3 ^ s4;
                next[1] = s3 ^ s5;
                next[2] = s0 ^ s4;
                next[3] = s0 ^ s1 ^ s5;
                next[4] = s1 ^ s2;
                next[5] = s2 ^ s3;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, p(x) = x^8 + x^4 + x^3 + x^2 + 1, multiply by alpha^8
                next[0] = s0 ^ s4 ^ s5 ^ s6;
                next[1] = s1 ^ s5 ^ s6 ^ s7;
                next[2] = s0 ^ s2 ^ s4 ^ s5 ^ s7;
                next[3] = s0 ^ s1 ^ s3 ^ s4;
                next[4] = s0 ^ s1 ^ s2 ^ s6;
                next[5] = s1 ^ s2 ^ s3 ^ s7;
                next[6] = s2 ^ s3 ^ s4;
                next[7] = s3 ^ s4 ^ s5;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, p(x) = x^10 + x^3 + 1, multiply by alpha^8
                next[0] = s2 ^ s9;
                next[1] = s3;
                next[2] = s4;
                next[3] = s2 ^ s5 ^ s9;
                next[4] = s3 ^ s6;
                next[5] = s4 ^ s7;
                next[6] = s5 ^ s8;
                next[7] = s6 ^ s9;
                next[8] = s0 ^ s7;
                next[9] = s1 ^ s8;
            end
        endcase
    end

    assign v_out = next;

endmodule

// Multiply a GF(2^m) element by alpha^16
// Supports m = 6, 8, 10
module chien_mul_a16 (
    input  [9:0] v_in,
    input  [1:0] m,
    output [9:0] v_out
);
    wire s0 = v_in[0];
    wire s1 = v_in[1];
    wire s2 = v_in[2];
    wire s3 = v_in[3];
    wire s4 = v_in[4];
    wire s5 = v_in[5];
    wire s6 = v_in[6];
    wire s7 = v_in[7];
    wire s8 = v_in[8];
    wire s9 = v_in[9];

    reg [9:0] next;

    always @* begin
        case (m)
            2'b01: begin
                // m = 6, multiply by alpha^16
                next[0] = s0 ^ s2 ^ s5;
                next[1] = s0 ^ s1 ^ s2 ^ s3 ^ s5;
                next[2] = s1 ^ s2 ^ s3 ^ s4;
                next[3] = s2 ^ s3 ^ s4 ^ s5;
                next[4] = s0 ^ s3 ^ s4 ^ s5;
                next[5] = s1 ^ s4 ^ s5;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, multiply by alpha^16
                next[0] = s2 ^ s5 ^ s7;
                next[1] = s3 ^ s6;
                next[2] = s0 ^ s2 ^ s4 ^ s5;
                next[3] = s0 ^ s1 ^ s2 ^ s3 ^ s6 ^ s7;
                next[4] = s1 ^ s3 ^ s4 ^ s5;
                next[5] = s2 ^ s4 ^ s5 ^ s6;
                next[6] = s0 ^ s3 ^ s5 ^ s6 ^ s7;
                next[7] = s1 ^ s4 ^ s6 ^ s7;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, multiply by alpha^16
                next[0] = s1 ^ s4 ^ s8;
                next[1] = s2 ^ s5 ^ s9;
                next[2] = s3 ^ s6;
                next[3] = s1 ^ s7 ^ s8;
                next[4] = s2 ^ s8 ^ s9;
                next[5] = s3 ^ s9;
                next[6] = s0 ^ s4;
                next[7] = s1 ^ s5;
                next[8] = s2 ^ s6;
                next[9] = s0 ^ s3 ^ s7;
            end
        endcase
    end

    assign v_out = next;

endmodule

// Multiply a GF(2^m) element by alpha^24
// Supports m = 6, 8, 10
module chien_mul_a24 (
    input  [9:0] v_in,
    input  [1:0] m,
    output [9:0] v_out
);
    wire s0 = v_in[0];
    wire s1 = v_in[1];
    wire s2 = v_in[2];
    wire s3 = v_in[3];
    wire s4 = v_in[4];
    wire s5 = v_in[5];
    wire s6 = v_in[6];
    wire s7 = v_in[7];
    wire s8 = v_in[8];
    wire s9 = v_in[9];

    reg [9:0] next;

    always @* begin
        case (m)
            2'b01: begin
                // m = 6, multiply by alpha^24
                next[0] = s0 ^ s2;
                next[1] = s1 ^ s2 ^ s3;
                next[2] = s2 ^ s3 ^ s4;
                next[3] = s3 ^ s4 ^ s5;
                next[4] = s0 ^ s4 ^ s5;
                next[5] = s1 ^ s5;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, multiply by alpha^24
                next[0] = s0 ^ s1;
                next[1] = s0 ^ s1 ^ s2;
                next[2] = s0 ^ s2 ^ s3;
                next[3] = s0 ^ s3 ^ s4;
                next[4] = s4 ^ s5;
                next[5] = s5 ^ s6;
                next[6] = s6 ^ s7;
                next[7] = s0 ^ s7;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, multiply by alpha^24
                next[0] = s0 ^ s6 ^ s7;
                next[1] = s1 ^ s7 ^ s8;
                next[2] = s2 ^ s8 ^ s9;
                next[3] = s0 ^ s3 ^ s6 ^ s7 ^ s9;
                next[4] = s0 ^ s1 ^ s4 ^ s7 ^ s8;
                next[5] = s1 ^ s2 ^ s5 ^ s8 ^ s9;
                next[6] = s2 ^ s3 ^ s6 ^ s9;
                next[7] = s3 ^ s4 ^ s7;
                next[8] = s4 ^ s5 ^ s8;
                next[9] = s5 ^ s6 ^ s9;
            end
        endcase
    end

    assign v_out = next;

endmodule


// Multiply a GF(2^m) element by alpha^32
// Supports m = 6, 8, 10
module chien_mul_a32 (
    input  [9:0] v_in,
    input  [1:0] m,
    output [9:0] v_out
);
    wire s0 = v_in[0];
    wire s1 = v_in[1];
    wire s2 = v_in[2];
    wire s3 = v_in[3];
    wire s4 = v_in[4];
    wire s5 = v_in[5];
    wire s6 = v_in[6];
    wire s7 = v_in[7];
    wire s8 = v_in[8];
    wire s9 = v_in[9];

    reg [9:0] next;

    always @* begin
        case (m)
            2'b01: begin
                // m = 6, multiply by alpha^32
                next[0] = s0 ^ s3;
                next[1] = s1 ^ s3 ^ s4;
                next[2] = s2 ^ s4 ^ s5;
                next[3] = s0 ^ s3 ^ s5;
                next[4] = s1 ^ s4;
                next[5] = s2 ^ s5;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, multiply by alpha^32
                next[0] = s0 ^ s1 ^ s4 ^ s7;
                next[1] = s1 ^ s2 ^ s5;
                next[2] = s0 ^ s1 ^ s2 ^ s3 ^ s4 ^ s6 ^ s7;
                next[3] = s0 ^ s2 ^ s3 ^ s5;
                next[4] = s0 ^ s3 ^ s6 ^ s7;
                next[5] = s1 ^ s4 ^ s7;
                next[6] = s2 ^ s5;
                next[7] = s0 ^ s3 ^ s6;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, multiply by alpha^32
                next[0] = s2 ^ s5 ^ s6 ^ s8;
                next[1] = s0 ^ s3 ^ s6 ^ s7 ^ s9;
                next[2] = s0 ^ s1 ^ s4 ^ s7 ^ s8;
                next[3] = s1 ^ s6 ^ s9;
                next[4] = s0 ^ s2 ^ s7;
                next[5] = s0 ^ s1 ^ s3 ^ s8;
                next[6] = s1 ^ s2 ^ s4 ^ s9;
                next[7] = s2 ^ s3 ^ s5;
                next[8] = s0 ^ s3 ^ s4 ^ s6;
                next[9] = s1 ^ s4 ^ s5 ^ s7;
            end
        endcase
    end

    assign v_out = next;

endmodule

module GF_mult(
	input [9:0] in_x1,
	input [9:0] in_x2,
	input [1:0] m,
	output [9:0] out_y
);

	//  mode 1: BCH(63,51) --> m = 6
	//  mode 2: BCH(255,239) --> m = 8
	//  mode 3: BCH(1023,983) --> m = 10

	// GF2 field multiplication: polynomial representation, XOR addition
	wire [18:0] mult;
	assign mult[0] = in_x1[0] & in_x2[0];
	assign mult[1] = (in_x1[1] & in_x2[0]) ^ (in_x1[0] & in_x2[1]);
	assign mult[2] = (in_x1[2] & in_x2[0]) ^ (in_x1[1] & in_x2[1]) ^ (in_x1[0] & in_x2[2]);
	assign mult[3] = (in_x1[3] & in_x2[0]) ^ (in_x1[2] & in_x2[1]) ^ (in_x1[1] & in_x2[2]) ^ (in_x1[0] & in_x2[3]);
	assign mult[4] = (in_x1[4] & in_x2[0]) ^ (in_x1[3] & in_x2[1]) ^ (in_x1[2] & in_x2[2]) ^ (in_x1[1] & in_x2[3]) ^ (in_x1[0] & in_x2[4]);
	assign mult[5] = (in_x1[5] & in_x2[0]) ^ (in_x1[4] & in_x2[1]) ^ (in_x1[3] & in_x2[2]) ^ (in_x1[2] & in_x2[3]) ^ (in_x1[1] & in_x2[4]) ^ (in_x1[0] & in_x2[5]);
	assign mult[6] = (in_x1[6] & in_x2[0]) ^ (in_x1[5] & in_x2[1]) ^ (in_x1[4] & in_x2[2]) ^ (in_x1[3] & in_x2[3]) ^ (in_x1[2] & in_x2[4]) ^ (in_x1[1] & in_x2[5]) ^ (in_x1[0] & in_x2[6]);
	assign mult[7] = (in_x1[7] & in_x2[0]) ^ (in_x1[6] & in_x2[1]) ^ (in_x1[5] & in_x2[2]) ^ (in_x1[4] & in_x2[3]) ^ (in_x1[3] & in_x2[4]) ^ (in_x1[2] & in_x2[5]) ^ (in_x1[1] & in_x2[6]) ^ (in_x1[0] & in_x2[7]);
	assign mult[8] = (in_x1[8] & in_x2[0]) ^ (in_x1[7] & in_x2[1]) ^ (in_x1[6] & in_x2[2]) ^ (in_x1[5] & in_x2[3]) ^ (in_x1[4] & in_x2[4]) ^ (in_x1[3] & in_x2[5]) ^ (in_x1[2] & in_x2[6]) ^ (in_x1[1] & in_x2[7]) ^ (in_x1[0] & in_x2[8]);
	assign mult[9] = (in_x1[9] & in_x2[0]) ^ (in_x1[8] & in_x2[1]) ^ (in_x1[7] & in_x2[2]) ^ (in_x1[6] & in_x2[3]) ^ (in_x1[5] & in_x2[4]) ^ (in_x1[4] & in_x2[5]) ^ (in_x1[3] & in_x2[6]) ^ (in_x1[2] & in_x2[7]) ^ (in_x1[1] & in_x2[8]) ^ (in_x1[0] & in_x2[9]);
	assign mult[10] = (in_x1[9] & in_x2[1]) ^ (in_x1[8] & in_x2[2]) ^ (in_x1[7] & in_x2[3]) ^ (in_x1[6] & in_x2[4]) ^ (in_x1[5] & in_x2[5]) ^ (in_x1[4] & in_x2[6]) ^ (in_x1[3] & in_x2[7]) ^ (in_x1[2] & in_x2[8]) ^ (in_x1[1] & in_x2[9]);
	assign mult[11] = (in_x1[9] & in_x2[2]) ^ (in_x1[8] & in_x2[3]) ^ (in_x1[7] & in_x2[4]) ^ (in_x1[6] & in_x2[5]) ^ (in_x1[5] & in_x2[6]) ^ (in_x1[4] & in_x2[7]) ^ (in_x1[3] & in_x2[8]) ^ (in_x1[2] & in_x2[9]);
	assign mult[12] = (in_x1[9] & in_x2[3]) ^ (in_x1[8] & in_x2[4]) ^ (in_x1[7] & in_x2[5]) ^ (in_x1[6] & in_x2[6]) ^ (in_x1[5] & in_x2[7]) ^ (in_x1[4] & in_x2[8]) ^ (in_x1[3] & in_x2[9]);
	assign mult[13] = (in_x1[9] & in_x2[4]) ^ (in_x1[8] & in_x2[5]) ^ (in_x1[7] & in_x2[6]) ^ (in_x1[6] & in_x2[7]) ^ (in_x1[5] & in_x2[8]) ^ (in_x1[4] & in_x2[9]);
	assign mult[14] = (in_x1[9] & in_x2[5]) ^ (in_x1[8] & in_x2[6]) ^ (in_x1[7] & in_x2[7]) ^ (in_x1[6] & in_x2[8]) ^ (in_x1[5] & in_x2[9]);
	assign mult[15] = (in_x1[9] & in_x2[6]) ^ (in_x1[8] & in_x2[7]) ^ (in_x1[7] & in_x2[8]) ^ (in_x1[6] & in_x2[9]);
	assign mult[16] = (in_x1[9] & in_x2[7]) ^ (in_x1[8] & in_x2[8]) ^ (in_x1[7] & in_x2[9]);
	assign mult[17] = (in_x1[9] & in_x2[8]) ^ (in_x1[8] & in_x2[9]);
	assign mult[18] = (in_x1[9] & in_x2[9]);

	reg [9:0] y;

	always @(*) begin
		case(m)
		2'b01: begin // reduction for m = 6;
			// reduction for m = 6; at most for x^5 * x^5 will have (5+5)= 10 degree
			// primitive polynomial: x^6 + x + 1 --> 10000011
			//x^6 = x^1 + x^0
			//x^7 = x^2 + x^1
			//x^8 = x^3 + x^2
			//x^9 = x^4 + x^3
			//x^10 = x^5 + x^4
			y[0] = mult[0] ^ mult[6];
			y[1] = mult[1] ^ mult[6] ^ mult[7];
			y[2] = mult[2] ^ mult[7] ^ mult[8];
			y[3] = mult[3] ^ mult[8] ^ mult[9];
			y[4] = mult[4] ^ mult[9] ^ mult[10];
			y[5] = mult[5] ^ mult[10];
			y[6] = 1'b0;
			y[7] = 1'b0;
			y[8] = 1'b0;
			y[9] = 1'b0;	
		end
		2'b10: begin // reduction for m = 8;
			// reduction for m = 8; at most will have x^7 * x^7 = x^14 degree
			// primitive polynomial: x^8 + x^4 + x^3 + x^2 + 1 --> 100011101
			//x^8 = x^4 + x^3 + x^2 + x^0
			//x^9 = x^5 + x^4 + x^3 + x^1
			//x^10 = x^6 + x^5 + x^4 + x^2
			//x^11 = x^7 + x^6 + x^5 + x^3
			//x^12 = x^8 + x^7 + x^6 + x^4 = x^4 + x^3 + x^2 + x^0 + x^7 + x^6 + x^4 = x^7 + x^6 + x^3 + x^2 + x^0
			//x^13 = x^9 + x^8 + x^7 + x^5 = x^5 + x^4 + x^3 + x^1 + x^4 + x^3 + x^2 + x^0 + x^7 + x^5 = x^7 + x^2 + x^1 + x^0
			//x^14 = x^13 * x = x^8 + x^3 + x^2 + x^1 = x^4 + x^3 + x^2 + x^0 + x^3 + x^2 + x^1 = x^4 + x^1 + x^0
			y[0] = mult[0] ^ mult[8] ^ mult[12] ^ mult[13] ^ mult[14];
			y[1] = mult[1] ^ mult[9] ^ mult[13] ^ mult[14];
			y[2] = mult[2] ^ mult[8] ^ mult[10] ^ mult[12] ^ mult[13];
			y[3] = mult[3] ^ mult[8] ^ mult[9] ^ mult[11] ^ mult[12] ;
			y[4] = mult[4] ^ mult[8] ^ mult[9] ^ mult[10] ^ mult[14];
			y[5] = mult[5] ^ mult[9] ^ mult[10] ^ mult[11];
			y[6] = mult[6] ^ mult[10] ^ mult[11] ^ mult[12];
			y[7] = mult[7] ^ mult[11] ^ mult[12] ^ mult[13];
			y[8] = 1'b0;	
			y[9] = 1'b0;
		end
		default: begin // reduction for m = 10;
			// reduction for m = 10; at most for x^9 * x^9 will have (9+9) = 18 degree
			// primitive polynomial: x^10 + x^3 + 1 --> 10000001001
			//x^10 = x^3 + x^0
			//x^11 = x^4 + x^1
			//x^12 = x^5 + x^2
			//x^13 = x^6 + x^3
			//x^14 = x^7 + x^4
			//x^15 = x^8 + x^5
			//x^16 = x^9 + x^6
			//x^17 = x^10 + x^7 = x^3 + x^0 + x^7
			//x^18 = x^11 + x^8 = x^4 + x^1 + x^8
			y[0] = mult[0] ^ mult[10] ^ mult[17];
			y[1] = mult[1] ^ mult[11] ^ mult[18];
			y[2] = mult[2] ^ mult[12];
			y[3] = mult[3] ^ mult[10] ^ mult[13] ^ mult[17];
			y[4] = mult[4] ^ mult[11] ^ mult[14] ^ mult[18];
			y[5] = mult[5] ^ mult[12] ^ mult[15];
			y[6] = mult[6] ^ mult[13] ^ mult[16];
			y[7] = mult[7] ^ mult[14] ^ mult[17];
			y[8] = mult[8] ^ mult[15] ^ mult[18];
			y[9] = mult[9] ^ mult[16];
		end
		endcase
	end	

	assign out_y = y;

endmodule

module MIN2_LLR(

    input clk,
    input rstn,
	input clear,

    input [6:0] abs_llr0,
	input [6:0] abs_llr1,
	input [6:0] abs_llr2,
	input [6:0] abs_llr3,
	input [6:0] abs_llr4,
	input [6:0] abs_llr5,
	input [6:0] abs_llr6,
	input [6:0] abs_llr7,
	input [6:0] i_min_val0,
	input [6:0] i_min_val1,
	input [9:0] i_min_idx0,
	input [9:0] i_min_idx1,

    input [7:0] counter,

    output [9:0] o_min_idx0,
    output [9:0] o_min_idx1,
	output [6:0] o_min_val0,
	output [6:0] o_min_val1
);

    wire [6:0] stage1_min_val [0:3];
    wire [9:0] stage1_min_idx [0:3];
    wire [6:0] stage1_max_val [0:3];
    wire [9:0] stage1_max_idx [0:3];

	wire [6:0] abs_llr [0:7];
	assign abs_llr[0] = abs_llr0;
	assign abs_llr[1] = abs_llr1;
	assign abs_llr[2] = abs_llr2;
	assign abs_llr[3] = abs_llr3;
	assign abs_llr[4] = abs_llr4;
	assign abs_llr[5] = abs_llr5;
	assign abs_llr[6] = abs_llr6;
	assign abs_llr[7] = abs_llr7;

	wire [9:0] index [0:7];
	assign index[0] = ((counter-1) << 3) + 0;
	assign index[1] = ((counter-1) << 3) + 1;
	assign index[2] = ((counter-1) << 3) + 2;
	assign index[3] = ((counter-1) << 3) + 3;
	assign index[4] = ((counter-1) << 3) + 4;
	assign index[5] = ((counter-1) << 3) + 5;
	assign index[6] = ((counter-1) << 3) + 6;
	assign index[7] = ((counter-1) << 3) + 7;
    genvar g;
    generate
        for (g = 0; g < 4; g = g + 1) begin : MINMAX1_BLOCK
            MINMAX2 MINMAX2_INST (
                .in_val0(abs_llr[2*g]),
                .in_idx0(index[2*g]),
                .in_val1(abs_llr[2*g + 1]),
                .in_idx1(index[2*g + 1]),
                .min_val(stage1_min_val[g]),
                .min_idx(stage1_min_idx[g]),
                .max_val(stage1_max_val[g]),
                .max_idx(stage1_max_idx[g])
            );
        end
    endgenerate

    wire [6:0] stage2_min_val [0:5];
    wire [9:0] stage2_min_idx [0:5];
    MINMAX2 MINMAX2_STAGE2_0 (
        .in_val0(stage1_min_val[0]),
        .in_idx0(stage1_min_idx[0]),
        .in_val1(stage1_min_val[1]),
        .in_idx1(stage1_min_idx[1]),
        .min_val(stage2_min_val[0]),
        .min_idx(stage2_min_idx[0]),
        .max_val(stage2_min_val[1]),
        .max_idx(stage2_min_idx[1])
    );
    MINMAX2 MINMAX2_STAGE2_1 (
        .in_val0(stage1_max_val[0]),
        .in_idx0(stage1_max_idx[0]),
        .in_val1(stage1_max_val[1]),
        .in_idx1(stage1_max_idx[1]),
        .min_val(stage2_min_val[2]),
        .min_idx(stage2_min_idx[2]),
        .max_val(),
        .max_idx()
    );
    MINMAX2 MINMAX2_STAGE2_3 (
        .in_val0(stage1_min_val[2]),
        .in_idx0(stage1_min_idx[2]),
        .in_val1(stage1_min_val[3]),
        .in_idx1(stage1_min_idx[3]),
        .min_val(stage2_min_val[3]),
        .min_idx(stage2_min_idx[3]),
        .max_val(stage2_min_val[4]),
        .max_idx(stage2_min_idx[4])
    );
    MINMAX2 MINMAX2_STAGE2_4 (
        .in_val0(stage1_max_val[2]),
        .in_idx0(stage1_max_idx[2]),
        .in_val1(stage1_max_val[3]),
        .in_idx1(stage1_max_idx[3]),
        .min_val(stage2_min_val[5]),
        .min_idx(stage2_min_idx[5]),
        .max_val(),
        .max_idx()
    );

    wire [6:0] stage3_min_val [0:3];
    wire [9:0] stage3_min_idx [0:3];

    reg [9:0] stage3_min_idx_w [0:3], stage3_min_idx_r [0:3];
    reg [6:0] stage3_min_val_w [0:3], stage3_min_val_r [0:3];
    always @(posedge clk) begin
		if(!rstn) begin
			stage3_min_idx_r[0] <= 10'd0;
			stage3_min_idx_r[1] <= 10'd0;
			stage3_min_idx_r[2] <= 10'd0;
			stage3_min_idx_r[3] <= 10'd0;
			stage3_min_val_r[0] <= 7'd0;
			stage3_min_val_r[1] <= 7'd0;
			stage3_min_val_r[2] <= 7'd0;
			stage3_min_val_r[3] <= 7'd0;
		end else begin
			stage3_min_idx_r[0] <= stage3_min_idx_w[0];
			stage3_min_idx_r[1] <= stage3_min_idx_w[1];
			stage3_min_idx_r[2] <= stage3_min_idx_w[2];
			stage3_min_idx_r[3] <= stage3_min_idx_w[3];
			stage3_min_val_r[0] <= stage3_min_val_w[0];
			stage3_min_val_r[1] <= stage3_min_val_w[1];
			stage3_min_val_r[2] <= stage3_min_val_w[2];
			stage3_min_val_r[3] <= stage3_min_val_w[3];
		end
    end

     always @(*) begin
		if(clear) begin
			stage3_min_idx_w[0] = 10'd0;
			stage3_min_idx_w[1] = 10'd0;
			stage3_min_idx_w[2] = 10'd0;
			stage3_min_idx_w[3] = 10'd0;
			stage3_min_val_w[0] = 7'd0;
			stage3_min_val_w[1] = 7'd0;
			stage3_min_val_w[2] = 7'd0;
			stage3_min_val_w[3] = 7'd0;
		end else  begin
			stage3_min_idx_w[0] = stage3_min_idx[0];
			stage3_min_idx_w[1] = stage3_min_idx[1];
			stage3_min_idx_w[2] = stage3_min_idx[2];
			stage3_min_idx_w[3] = stage3_min_idx[3];
			stage3_min_val_w[0] = stage3_min_val[0];
			stage3_min_val_w[1] = stage3_min_val[1];
			stage3_min_val_w[2] = stage3_min_val[2];
			stage3_min_val_w[3] = stage3_min_val[3];
		end
    end
    MINMAX2 MINMAX2_STAGE3_0 (
        .in_val0(stage2_min_val[0]),
        .in_idx0(stage2_min_idx[0]),
        .in_val1(stage2_min_val[3]),
        .in_idx1(stage2_min_idx[3]),
        .min_val(stage3_min_val[0]),
        .min_idx(stage3_min_idx[0]),
        .max_val(stage3_min_val[1]),
        .max_idx(stage3_min_idx[1])
    );
    MINMAX2 MINMAX2_STAGE3_1 (
        .in_val0(stage2_min_val[1]),
        .in_idx0(stage2_min_idx[1]),
        .in_val1(stage2_min_val[2]),
        .in_idx1(stage2_min_idx[2]),
        .min_val(stage3_min_val[2]),
        .min_idx(stage3_min_idx[2]),
        .max_val(),
        .max_idx()
    );
    MINMAX2 MINMAX2_STAGE3_2 (
        .in_val0(stage2_min_val[4]),
        .in_idx0(stage2_min_idx[4]),
        .in_val1(stage2_min_val[5]),
        .in_idx1(stage2_min_idx[5]),
        .min_val(stage3_min_val[3]),
        .min_idx(stage3_min_idx[3]),
        .max_val(),
        .max_idx()
    );

    wire [6:0] stage4_min_val [0:2];
    wire [9:0] stage4_min_idx [0:2];
    MINMAX2 MINMAX2_STAGE4_0 (
        .in_val0(i_min_val0),
        .in_idx0(i_min_idx0),
        .in_val1(stage3_min_val_r[0]),
        .in_idx1(stage3_min_idx_r[0]),
        .min_val(o_min_val0),
        .min_idx(o_min_idx0),
        .max_val(stage4_min_val[0]),
        .max_idx(stage4_min_idx[0])
    );
    MINMAX2 MINMAX2_STAGE4_1 (
        .in_val0(stage3_min_val_r[1]),
        .in_idx0(stage3_min_idx_r[1]),
        .in_val1(stage3_min_val_r[2]),
        .in_idx1(stage3_min_idx_r[2]),
        .min_val(stage4_min_val[1]),
        .min_idx(stage4_min_idx[1]),
        .max_val(),
        .max_idx()
    );
    MINMAX2 MINMAX2_STAGE4_2 (
        .in_val0(stage3_min_val_r[3]),
        .in_idx0(stage3_min_idx_r[3]),
        .in_val1(i_min_val1),
        .in_idx1(i_min_idx1),
        .min_val(stage4_min_val[2]),
        .min_idx(stage4_min_idx[2]),
        .max_val(),
        .max_idx()
    );

    wire [6:0] stage5_min_val;
    wire [9:0] stage5_min_idx;
    MINMAX2 MINMAX2_STAGE5 (
        .in_val0(stage4_min_val[0]),
        .in_idx0(stage4_min_idx[0]),
        .in_val1(stage4_min_val[1]),
        .in_idx1(stage4_min_idx[1]),
        .min_val(stage5_min_val),
        .min_idx(stage5_min_idx),
        .max_val(),
        .max_idx()
    );

    MINMAX2 MINMAX2_STAGE6 (
        .in_val0(stage5_min_val),
        .in_idx0(stage5_min_idx),
        .in_val1(stage4_min_val[2]),
        .in_idx1(stage4_min_idx[2]),
        .min_val(o_min_val1),
        .min_idx(o_min_idx1),
        .max_val(),
        .max_idx()
    );
endmodule

module MINMAX2(
    input [6:0] in_val0,
    input [9:0] in_idx0,
    input [6:0] in_val1,
    input [9:0] in_idx1,

    output [6:0] min_val,
    output [9:0] min_idx,
    output [6:0] max_val,
    output [9:0] max_idx
);
    assign {min_val, min_idx} = (in_val0 <= in_val1) ? {in_val0, in_idx0} : {in_val1, in_idx1};
    assign {max_val, max_idx} = (in_val0 > in_val1) ? {in_val0, in_idx0} : {in_val1, in_idx1};
endmodule


module SYNDROME_GEN(
	input [9:0] i_syndrome_1,
	input [9:0] i_syndrome_2,
	input [9:0] i_syndrome_3,
	input [9:0] i_syndrome_4,
	input [9:0] i_syndrome_5,
	input [9:0] i_syndrome_6,
	input [9:0] i_syndrome_7,
	input [9:0] i_syndrome_8,
	input [7:0] syndrome_data_in,
	input [1:0] code,
	output [9:0] o_syndrome_1,
	output [9:0] o_syndrome_2,
	output [9:0] o_syndrome_3,
	output [9:0] o_syndrome_4,
	output [9:0] o_syndrome_5,
	output [9:0] o_syndrome_6,
	output [9:0] o_syndrome_7,
	output [9:0] o_syndrome_8
);

wire [9:0] syndrome [0:7];
assign syndrome[0] = i_syndrome_1;
assign syndrome[1] = i_syndrome_2;
assign syndrome[2] = i_syndrome_3;
assign syndrome[3] = i_syndrome_4;
assign syndrome[4] = i_syndrome_5;
assign syndrome[5] = i_syndrome_6;
assign syndrome[6] = i_syndrome_7;
assign syndrome[7] = i_syndrome_8;

wire [9:0] syndrome_out [0:7];
assign o_syndrome_1 = syndrome_out[0];
assign o_syndrome_2 = syndrome_out[1];
assign o_syndrome_3 = syndrome_out[2];
assign o_syndrome_4 = syndrome_out[3];
assign o_syndrome_5 = syndrome_out[4];
assign o_syndrome_6 = syndrome_out[5];
assign o_syndrome_7 = syndrome_out[6];
assign o_syndrome_8 = syndrome_out[7];

syndrome_block_n1 SYNDROME_BLOCK_N1_INST (
	.S_in(syndrome[0]),
	.x(syndrome_data_in),
	.m(code),
	.S_out(syndrome_out[0])
);
syndrome_block_n2 SYNDROME_BLOCK_N2_INST (
	.S_in(syndrome[1]),
	.x(syndrome_data_in),
	.m(code),
	.S_out(syndrome_out[1])
);
syndrome_block_n3 SYNDROME_BLOCK_N3_INST (
	.S_in(syndrome[2]),
	.x(syndrome_data_in),
	.m(code),
	.S_out(syndrome_out[2])
);
syndrome_block_n4 SYNDROME_BLOCK_N4_INST (
	.S_in(syndrome[3]),
	.x(syndrome_data_in),
	.m(code),
	.S_out(syndrome_out[3])
);
syndrome_block_n5 SYNDROME_BLOCK_N5_INST (
	.S_in(syndrome[4]),
	.x(syndrome_data_in),
	.S_out(syndrome_out[4])
);
syndrome_block_n6 SYNDROME_BLOCK_N6_INST (
	.S_in(syndrome[5]),
	.x(syndrome_data_in),
	.S_out(syndrome_out[5])
);
syndrome_block_n7 SYNDROME_BLOCK_N7_INST (
	.S_in(syndrome[6]),
	.x(syndrome_data_in),
	.S_out(syndrome_out[6])
);
syndrome_block_n8 SYNDROME_BLOCK_N8_INST (
	.S_in(syndrome[7]),	
	.x(syndrome_data_in),
	.S_out(syndrome_out[7])
);


endmodule

// Syndrome block update for n = 1 (S1)
// Supports m = 6, 8, 10 via case(m)
module syndrome_block_n1 (
    input  [9:0] S_in,   // current S_n(i)
    input  [7:0] x,      // 8 input bits, x[0]=x0 ... x[7]=x7
    input  [1:0] m,      // 6, 8, or 10
    output [9:0] S_out   // next S_n(i+8)
);
    // Alias bits for readability
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        next = 10'b0;
        case (m)
            2'b01: begin
                // m = 6, n = 1
                next[0] = s3 ^ s4 ^ x1 ^ x7;
                next[1] = s3 ^ s5 ^ x0 ^ x1 ^ x6;
                next[2] = s0 ^ s4 ^ x0 ^ x5;
                next[3] = s0 ^ s1 ^ s5 ^ x4;
                next[4] = s1 ^ s2 ^ x3;
                next[5] = s2 ^ s3 ^ x2;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, n = 1
                next[0] = s0 ^ s4 ^ s5 ^ s6 ^ x7;
                next[1] = s1 ^ s5 ^ s6 ^ s7 ^ x6;
                next[2] = s0 ^ s2 ^ s4 ^ s5 ^ s7 ^ x5;
                next[3] = s0 ^ s1 ^ s3 ^ s4 ^ x4;
                next[4] = s0 ^ s1 ^ s2 ^ s6 ^ x3;
                next[5] = s1 ^ s2 ^ s3 ^ s7 ^ x2;
                next[6] = s2 ^ s3 ^ s4 ^ x1;
                next[7] = s3 ^ s4 ^ s5 ^ x0;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, n = 1
                next[0] = s2 ^ s9 ^ x7;
                next[1] = s3 ^ x6;
                next[2] = s4 ^ x5;
                next[3] = s2 ^ s5 ^ s9 ^ x4;
                next[4] = s3 ^ s6 ^ x3;
                next[5] = s4 ^ s7 ^ x2;
                next[6] = s5 ^ s8 ^ x1;
                next[7] = s6 ^ s9 ^ x0;
                next[8] = s0 ^ s7;
                next[9] = s1 ^ s8;
            end
        endcase
    end

    assign S_out = next;
endmodule

// Syndrome block update for n = 2 (S2)
module syndrome_block_n2 (
    input  [9:0] S_in,
    input  [7:0] x,
    input  [1:0] m,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        next = 10'b0;
        case (m)
            2'b01: begin
                // m = 6, n = 2
                next[0] = s0 ^ s2 ^ s5 ^ x1 ^ x4 ^ x7;
                next[1] = s0 ^ s1 ^ s2 ^ s3 ^ s5 ^ x4;
                next[2] = s1 ^ s2 ^ s3 ^ s4 ^ x0 ^ x1 ^ x3 ^ x6;
                next[3] = s2 ^ s3 ^ s4 ^ s5 ^ x3;
                next[4] = s0 ^ s3 ^ s4 ^ s5 ^ x0 ^ x2 ^ x5;
                next[5] = s1 ^ s4 ^ s5 ^ x2;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, n = 2
                next[0] = s2 ^ s5 ^ s7 ^ x0 ^ x1 ^ x3 ^ x7;
                next[1] = s3 ^ s6 ^ x0;
                next[2] = s0 ^ s2 ^ s4 ^ s5 ^ x1 ^ x2 ^ x3 ^ x6;
                next[3] = s0 ^ s1 ^ s2 ^ s3 ^ s6 ^ s7 ^ x1 ^ x3;
                next[4] = s1 ^ s3 ^ s4 ^ s5 ^ x0 ^ x2 ^ x3 ^ x5;
                next[5] = s2 ^ s4 ^ s5 ^ s6 ^ x2;
                next[6] = s0 ^ s3 ^ s5 ^ s6 ^ s7 ^ x1 ^ x2 ^ x4;
                next[7] = s1 ^ s4 ^ s6 ^ s7 ^ x1;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, n = 2
                next[0] = s1 ^ s4 ^ s8 ^ x2 ^ x7;
                next[1] = s2 ^ s5 ^ s9;
                next[2] = s3 ^ s6 ^ x1 ^ x6;
                next[3] = s1 ^ s7 ^ s8 ^ x2;
                next[4] = s2 ^ s8 ^ s9 ^ x0 ^ x5;
                next[5] = s3 ^ s9 ^ x1;
                next[6] = s0 ^ s4 ^ x4;
                next[7] = s1 ^ s5 ^ x0;
                next[8] = s2 ^ s6 ^ x3;
                next[9] = s0 ^ s3 ^ s7;
            end
        endcase
    end

    assign S_out = next;

endmodule
// Syndrome block update for n = 3 (S3)
module syndrome_block_n3 (
    input  [9:0] S_in,
    input  [7:0] x,
    input  [1:0] m,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];
    reg [9:0] next;

    always @* begin
        next = 10'b0;
        case (m)
            2'b01: begin
                // m = 6, n = 3
                next[0] = s0 ^ s2 ^ x0 ^ x1 ^ x3 ^ x5 ^ x7;
                next[1] = s1 ^ s2 ^ s3 ^ x0 ^ x1 ^ x5;
                next[2] = s2 ^ s3 ^ s4 ^ x1 ^ x3;
                next[3] = s3 ^ s4 ^ s5 ^ x0 ^ x1 ^ x2 ^ x4 ^ x6;
                next[4] = s0 ^ s4 ^ s5 ^ x0 ^ x4;
                next[5] = s1 ^ s5 ^ x0 ^ x2;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, n = 3
                next[0] = s0 ^ s1 ^ x0 ^ x1 ^ x3 ^ x7;
                next[1] = s0 ^ s1 ^ s2 ^ x2 ^ x4;
                next[2] = s0 ^ s2 ^ s3 ^ x0 ^ x1 ^ x2 ^ x3;
                next[3] = s0 ^ s3 ^ s4 ^ x1 ^ x3 ^ x4 ^ x6;
                next[4] = s4 ^ s5 ^ x0 ^ x4;
                next[5] = s5 ^ s6 ^ x0 ^ x1 ^ x2 ^ x4;
                next[6] = s6 ^ s7 ^ x0 ^ x3 ^ x5;
                next[7] = s0 ^ s7 ^ x3;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, n = 3
                next[0] = s0 ^ s6 ^ s7 ^ x7;
                next[1] = s1 ^ s7 ^ s8 ^ x0 ^ x1;
                next[2] = s2 ^ s8 ^ s9 ^ x3;
                next[3] = s0 ^ s3 ^ s6 ^ s7 ^ s9 ^ x6;
                next[4] = s0 ^ s1 ^ s4 ^ s7 ^ s8 ^ x1;
                next[5] = s1 ^ s2 ^ s5 ^ s8 ^ s9 ^ x2 ^ x3;
                next[6] = s2 ^ s3 ^ s6 ^ s9 ^ x5;
                next[7] = s3 ^ s4 ^ s7 ^ x0;
                next[8] = s4 ^ s5 ^ s8 ^ x1 ^ x2;
                next[9] = s5 ^ s6 ^ s9 ^ x4;
            end
        endcase
    end

    assign S_out = next;

endmodule
// Syndrome block update for n = 4 (S4)
module syndrome_block_n4 (
    input  [9:0] S_in,
    input  [7:0] x,
    input  [1:0] m,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        next = 10'b0;
        case (m)
            2'b01: begin
                // m = 6, n = 4
                next[0] = s0 ^ s3 ^ x1 ^ x3 ^ x4 ^ x7;
                next[1] = s1 ^ s3 ^ s4 ^ x3;
                next[2] = s2 ^ s4 ^ s5 ^ x0 ^ x2 ^ x4 ^ x5;
                next[3] = s0 ^ s3 ^ s5 ^ x0 ^ x2 ^ x5;
                next[4] = s1 ^ s4 ^ x0 ^ x1 ^ x2 ^ x3 ^ x6;
                next[5] = s2 ^ s5 ^ x2;
                next[9:6] = 4'b0;
            end
            2'b10: begin
                // m = 8, n = 4
                next[0] = s0 ^ s1 ^ s4 ^ s7 ^ x1 ^ x4 ^ x5 ^ x7;
                next[1] = s1 ^ s2 ^ s5 ^ x1;
                next[2] = s0 ^ s1 ^ s2 ^ s3 ^ s4 ^ s6 ^ s7 ^
                          x1 ^ x2 ^ x3 ^ x4 ^ x5;
                next[3] = s0 ^ s2 ^ s3 ^ s5 ^ x0 ^ x1 ^ x3 ^ x4 ^ x5;
                next[4] = s0 ^ s3 ^ s6 ^ s7 ^ x0 ^ x2 ^ x5 ^ x6;
                next[5] = s1 ^ s4 ^ s7 ^ x2;
                next[6] = s2 ^ s5 ^ x3 ^ x4;
                next[7] = s0 ^ s3 ^ s6 ^ x1 ^ x2 ^ x4;
                next[9:8] = 2'b0;
            end
            default: begin
                // m = 10, n = 4
                next[0] = s2 ^ s5 ^ s6 ^ s8 ^ x1 ^ x2 ^ x7;
                next[1] = s0 ^ s3 ^ s6 ^ s7 ^ s9;
                next[2] = s0 ^ s1 ^ s4 ^ s7 ^ s8 ^ x4;
                next[3] = s1 ^ s6 ^ s9 ^ x1;
                next[4] = s0 ^ s2 ^ s7 ^ x0 ^ x1 ^ x6;
                next[5] = s0 ^ s1 ^ s3 ^ s8 ^ x4;
                next[6] = s1 ^ s2 ^ s4 ^ s9 ^ x2 ^ x3;
                next[7] = s2 ^ s3 ^ s5 ^ x0;
                next[8] = s0 ^ s3 ^ s4 ^ s6 ^ x0 ^ x5;
                next[9] = s1 ^ s4 ^ s5 ^ s7 ^ x3;
            end
        endcase
    end

    assign S_out = next;

endmodule

// Syndrome block update for n = 5, m = 10 only
module syndrome_block_n5 (
    input  [9:0] S_in,
    input  [7:0] x,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        // m = 10, n = 5
        next[0] = s0 ^ s5 ^ s8 ^ x1 ^ x3 ^ x5 ^ x7;
        next[1] = s1 ^ s6 ^ s9 ^ x0 ^ x2;
        next[2] = s0 ^ s2 ^ s7;
        next[3] = s1 ^ s3 ^ s5 ^ x1 ^ x5;
        next[4] = s2 ^ s4 ^ s6 ^ x2;
        next[5] = s0 ^ s3 ^ s5 ^ s7 ^ x0 ^ x2 ^ x4 ^ x6;
        next[6] = s1 ^ s4 ^ s6 ^ s8 ^ x1 ^ x3;
        next[7] = s2 ^ s5 ^ s7 ^ s9 ^ x0;
        next[8] = s3 ^ s6 ^ s8 ^ x0 ^ x4;
        next[9] = s4 ^ s7 ^ s9 ^ x1;
    end

    assign S_out = next;

endmodule

// Syndrome block update for n = 6, m = 10 only
module syndrome_block_n6 (
    input  [9:0] S_in,
    input  [7:0] x,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        // m = 10, n = 6
        next[0] = s0 ^ s2 ^ s4 ^ s9 ^ x2 ^ x3 ^ x7;
        next[1] = s1 ^ s3 ^ s5 ^ x4;
        next[2] = s2 ^ s4 ^ s6 ^ x0 ^ x1 ^ x5;
        next[3] = s2 ^ s3 ^ s4 ^ s5 ^ s7 ^ s9 ^ x2 ^ x3;
        next[4] = s3 ^ s4 ^ s5 ^ s6 ^ s8 ^ x0 ^ x3 ^ x4;
        next[5] = s4 ^ s5 ^ s6 ^ s7 ^ s9 ^ x5;
        next[6] = s0 ^ s5 ^ s6 ^ s7 ^ s8 ^ x1 ^ x2 ^ x6;
        next[7] = s1 ^ s6 ^ s7 ^ s8 ^ s9 ^ x0;
        next[8] = s0 ^ s2 ^ s7 ^ s8 ^ s9 ^ x1 ^ x4;
        next[9] = s1 ^ s3 ^ s8 ^ s9 ^ x1 ^ x2;
    end

    assign S_out = next;

endmodule

// Syndrome block update for n = 7, m = 10 only
module syndrome_block_n7 (
    input  [9:0] S_in,
    input  [7:0] x,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        // m = 10, n = 7
        next[0] = s1 ^ s2 ^ s3 ^ s4 ^ s6 ^ s8 ^ s9 ^ x7;
        next[1] = s2 ^ s3 ^ s4 ^ s5 ^ s7 ^ s9 ^ x0 ^ x2 ^ x4;
        next[2] = s3 ^ s4 ^ s5 ^ s6 ^ s8 ^ x1;
        next[3] = s1 ^ s2 ^ s3 ^ s5 ^ s7 ^ s8;
        next[4] = s0 ^ s2 ^ s3 ^ s4 ^ s6 ^ s8 ^ s9 ^ x1 ^ x3 ^ x5;
        next[5] = s1 ^ s3 ^ s4 ^ s5 ^ s7 ^ s9 ^ x2;
        next[6] = s0 ^ s2 ^ s4 ^ s5 ^ s6 ^ s8;
        next[7] = s0 ^ s1 ^ s3 ^ s5 ^ s6 ^ s7 ^ s9 ^
                  x0 ^ x1 ^ x2 ^ x3 ^ x4 ^ x5 ^ x6;
        next[8] = s0 ^ s1 ^ s2 ^ s4 ^ s6 ^ s7 ^ s8 ^ x2 ^ x3;
        next[9] = s0 ^ s1 ^ s2 ^ s3 ^ s5 ^ s7 ^ s8 ^ s9 ^ x0;
    end

    assign S_out = next;

endmodule

// Syndrome block update for n = 8, m = 10 only
module syndrome_block_n8 (
    input  [9:0] S_in,
    input  [7:0] x,
    output [9:0] S_out
);
    wire s0 = S_in[0];
    wire s1 = S_in[1];
    wire s2 = S_in[2];
    wire s3 = S_in[3];
    wire s4 = S_in[4];
    wire s5 = S_in[5];
    wire s6 = S_in[6];
    wire s7 = S_in[7];
    wire s8 = S_in[8];
    wire s9 = S_in[9];

    wire x0 = x[7];
    wire x1 = x[6];
    wire x2 = x[5];
    wire x3 = x[4];
    wire x4 = x[3];
    wire x5 = x[2];
    wire x6 = x[1];
    wire x7 = x[0];

    reg [9:0] next;

    always @* begin
        // m = 10, n = 8
        next[0] = s0 ^ s1 ^ s2 ^ s4 ^ s6 ^ s7 ^ s9 ^
                  x1 ^ x2 ^ x4 ^ x7;
        next[1] = s1 ^ s2 ^ s3 ^ s5 ^ s7 ^ s8 ^ x3;
        next[2] = s0 ^ s2 ^ s3 ^ s4 ^ s6 ^ s8 ^ s9 ^ x2 ^ x3;
        next[3] = s0 ^ s2 ^ s3 ^ s5 ^ s6 ^ x4;
        next[4] = s0 ^ s1 ^ s3 ^ s4 ^ s6 ^ s7 ^
                  x0 ^ x3 ^ x4;
        next[5] = s1 ^ s2 ^ s4 ^ s5 ^ s7 ^ s8 ^ x2 ^ x3;
        next[6] = s0 ^ s2 ^ s3 ^ s5 ^ s6 ^ s8 ^ s9 ^
                  x0 ^ x1 ^ x5;
        next[7] = s1 ^ s3 ^ s4 ^ s6 ^ s7 ^ s9 ^ x0;
        next[8] = s0 ^ s2 ^ s4 ^ s5 ^ s7 ^ s8 ^
                  x0 ^ x1 ^ x3 ^ x6;
        next[9] = s0 ^ s1 ^ s3 ^ s5 ^ s6 ^ s8 ^ s9 ^
                  x0 ^ x5;
    end

    assign S_out = next;

endmodule


module BERLEKAMP_MULT(
	input [9:0] i_gf_mult_x1_0,
	input [9:0] i_gf_mult_x1_1,
	input [9:0] i_gf_mult_x1_2,
	input [9:0] i_gf_mult_x1_3,
	input [9:0] i_gf_mult_x1_4,
	input [9:0] i_gf_mult_x1_5,
	input [9:0] i_gf_mult_x1_6,
	input [9:0] i_gf_mult_x1_7,
	input [9:0] i_gf_mult_x1_8,
	input [9:0] i_gf_mult_x1_9,
	input [9:0] i_gf_mult_x2_0,
	input [9:0] i_gf_mult_x2_1,
	input [9:0] i_gf_mult_x2_2,
	input [9:0] i_gf_mult_x2_3,
	input [9:0] i_gf_mult_x2_4,
	input [9:0] i_gf_mult_x2_5,
	input [9:0] i_gf_mult_x2_6,
	input [9:0] i_gf_mult_x2_7,
	input [9:0] i_gf_mult_x2_8,
	input [9:0] i_gf_mult_x2_9,
	output [9:0] o_gf_mult_y0,
	output [9:0] o_gf_mult_y1,
	output [9:0] o_gf_mult_y2,
	output [9:0] o_gf_mult_y3,
	output [9:0] o_gf_mult_y4,
	output [9:0] o_gf_mult_y5,
	output [9:0] o_gf_mult_y6,
	output [9:0] o_gf_mult_y7,
	output [9:0] o_gf_mult_y8,
	output [9:0] o_gf_mult_y9,
	input [1:0] code
);

wire [9:0] gf_mult_in1 [0:9];
wire [9:0] gf_mult_in2 [0:9];
wire [9:0] gf_mult_out [0:9];

assign gf_mult_in1[0] = i_gf_mult_x1_0;
assign gf_mult_in1[1] = i_gf_mult_x1_1;
assign gf_mult_in1[2] = i_gf_mult_x1_2;
assign gf_mult_in1[3] = i_gf_mult_x1_3;
assign gf_mult_in1[4] = i_gf_mult_x1_4;
assign gf_mult_in1[5] = i_gf_mult_x1_5;
assign gf_mult_in1[6] = i_gf_mult_x1_6;
assign gf_mult_in1[7] = i_gf_mult_x1_7;
assign gf_mult_in1[8] = i_gf_mult_x1_8;
assign gf_mult_in1[9] = i_gf_mult_x1_9;

assign gf_mult_in2[0] = i_gf_mult_x2_0;
assign gf_mult_in2[1] = i_gf_mult_x2_1;
assign gf_mult_in2[2] = i_gf_mult_x2_2;
assign gf_mult_in2[3] = i_gf_mult_x2_3;
assign gf_mult_in2[4] = i_gf_mult_x2_4;
assign gf_mult_in2[5] = i_gf_mult_x2_5	;
assign gf_mult_in2[6] = i_gf_mult_x2_6;
assign gf_mult_in2[7] = i_gf_mult_x2_7;
assign gf_mult_in2[8] = i_gf_mult_x2_8;
assign gf_mult_in2[9] = i_gf_mult_x2_9;

assign o_gf_mult_y0 = gf_mult_out[0];
assign o_gf_mult_y1 = gf_mult_out[1];
assign o_gf_mult_y2 = gf_mult_out[2];
assign o_gf_mult_y3 = gf_mult_out[3];
assign o_gf_mult_y4 = gf_mult_out[4];
assign o_gf_mult_y5 = gf_mult_out[5];
assign o_gf_mult_y6 = gf_mult_out[6];
assign o_gf_mult_y7 = gf_mult_out[7];
assign o_gf_mult_y8 = gf_mult_out[8];
assign o_gf_mult_y9 = gf_mult_out[9];
genvar g;
generate
	for (g = 0; g < 10; g = g + 1) begin : GF_MULT_BLOCK
		GF_mult GF_MULT_INST (
			.in_x1(gf_mult_in1[g]),
			.in_x2(gf_mult_in2[g]),
			.m(code),
			.out_y(gf_mult_out[g])
		);
	end
endgenerate

endmodule

module MINCORR(
    input [9:0] in_corr0,
    input [9:0] in_corr1,
    input [9:0] in_corr2,
    input [9:0] in_corr3,

    output [1:0] min_corr
);
    wire [1:0] min_idx0;
    wire [9:0] min_val0;
    MIN2 MIN2_INST1(
        .in_val0(in_corr0),
        .in_idx0(2'b00),
        .in_val1(in_corr1),
        .in_idx1(2'b01),
        .min_val(min_val0),
        .min_idx(min_idx0)
    );
    wire [1:0] min_idx1;
    wire [9:0] min_val1;
    MIN2 MIN2_INST2(
        .in_val0(in_corr2),
        .in_idx0(2'b10),
        .in_val1(in_corr3),
        .in_idx1(2'b11),
        .min_val(min_val1),
        .min_idx(min_idx1)
    );
    MIN2 MIN2_INST3(
        .in_val0(min_val0),
        .in_idx0(min_idx0),
        .in_val1(min_val1),
        .in_idx1(min_idx1),
        .min_val(),
        .min_idx(min_corr)
    );
endmodule

module MIN2(
    input [9:0] in_val0,
    input [1:0] in_idx0,
    input [9:0] in_val1,
    input [1:0] in_idx1,

    output [9:0] min_val,
    output [1:0] min_idx
);
    assign {min_val, min_idx} = (in_val0 <= in_val1) ? {in_val0, in_idx0} : {in_val1, in_idx1};
endmodule



// Auto-generated GF(2^m) multiply-by-alpha^k modules
// m encoding:
//   2'b01 -> m = 6  (GF(2^6))
//   2'b10 -> m = 8  (GF(2^8))
//   2'b11 -> m = 10 (GF(2^10))
//   default -> treated as m = 8 (same as 2'b10)

// ==============================================
// GF multiply by alpha^1
// ==============================================
module gf_mul_alpha_1 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[5];
                next[1] = current[0] ^ current[5];
                next[2] = current[1];
                next[3] = current[2];
                next[4] = current[3];
                next[5] = current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[7];
                next[1] = current[0];
                next[2] = current[1] ^ current[7];
                next[3] = current[2] ^ current[7];
                next[4] = current[3] ^ current[7];
                next[5] = current[4];
                next[6] = current[5];
                next[7] = current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[9];
                next[1] = current[0];
                next[2] = current[1];
                next[3] = current[2] ^ current[9];
                next[4] = current[3];
                next[5] = current[4];
                next[6] = current[5];
                next[7] = current[6];
                next[8] = current[7];
                next[9] = current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[7];
                next[1] = current[0];
                next[2] = current[1] ^ current[7];
                next[3] = current[2] ^ current[7];
                next[4] = current[3] ^ current[7];
                next[5] = current[4];
                next[6] = current[5];
                next[7] = current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^2
// ==============================================
module gf_mul_alpha_2 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[4];
                next[1] = current[4] ^ current[5];
                next[2] = current[0] ^ current[5];
                next[3] = current[1];
                next[4] = current[2];
                next[5] = current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[6];
                next[1] = current[7];
                next[2] = current[0] ^ current[6];
                next[3] = current[1] ^ current[6] ^ current[7];
                next[4] = current[2] ^ current[6] ^ current[7];
                next[5] = current[3] ^ current[7];
                next[6] = current[4];
                next[7] = current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[8];
                next[1] = current[9];
                next[2] = current[0];
                next[3] = current[1] ^ current[8];
                next[4] = current[2] ^ current[9];
                next[5] = current[3];
                next[6] = current[4];
                next[7] = current[5];
                next[8] = current[6];
                next[9] = current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[6];
                next[1] = current[7];
                next[2] = current[0] ^ current[6];
                next[3] = current[1] ^ current[6] ^ current[7];
                next[4] = current[2] ^ current[6] ^ current[7];
                next[5] = current[3] ^ current[7];
                next[6] = current[4];
                next[7] = current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^3
// ==============================================
module gf_mul_alpha_3 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[3];
                next[1] = current[3] ^ current[4];
                next[2] = current[4] ^ current[5];
                next[3] = current[0] ^ current[5];
                next[4] = current[1];
                next[5] = current[2];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[5];
                next[1] = current[6];
                next[2] = current[5] ^ current[7];
                next[3] = current[0] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[5] = current[2] ^ current[6] ^ current[7];
                next[6] = current[3] ^ current[7];
                next[7] = current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[7];
                next[1] = current[8];
                next[2] = current[9];
                next[3] = current[0] ^ current[7];
                next[4] = current[1] ^ current[8];
                next[5] = current[2] ^ current[9];
                next[6] = current[3];
                next[7] = current[4];
                next[8] = current[5];
                next[9] = current[6];
            end
            default: begin // default to GF(2^8)
                next[0] = current[5];
                next[1] = current[6];
                next[2] = current[5] ^ current[7];
                next[3] = current[0] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[5] = current[2] ^ current[6] ^ current[7];
                next[6] = current[3] ^ current[7];
                next[7] = current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^4
// ==============================================
module gf_mul_alpha_4 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[2];
                next[1] = current[2] ^ current[3];
                next[2] = current[3] ^ current[4];
                next[3] = current[4] ^ current[5];
                next[4] = current[0] ^ current[5];
                next[5] = current[1];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[4];
                next[1] = current[5];
                next[2] = current[4] ^ current[6];
                next[3] = current[4] ^ current[5] ^ current[7];
                next[4] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[5] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[6] = current[2] ^ current[6] ^ current[7];
                next[7] = current[3] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[6];
                next[1] = current[7];
                next[2] = current[8];
                next[3] = current[6] ^ current[9];
                next[4] = current[0] ^ current[7];
                next[5] = current[1] ^ current[8];
                next[6] = current[2] ^ current[9];
                next[7] = current[3];
                next[8] = current[4];
                next[9] = current[5];
            end
            default: begin // default to GF(2^8)
                next[0] = current[4];
                next[1] = current[5];
                next[2] = current[4] ^ current[6];
                next[3] = current[4] ^ current[5] ^ current[7];
                next[4] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[5] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[6] = current[2] ^ current[6] ^ current[7];
                next[7] = current[3] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^5
// ==============================================
module gf_mul_alpha_5 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1];
                next[1] = current[1] ^ current[2];
                next[2] = current[2] ^ current[3];
                next[3] = current[3] ^ current[4];
                next[4] = current[4] ^ current[5];
                next[5] = current[0] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[3] ^ current[7];
                next[1] = current[4];
                next[2] = current[3] ^ current[5] ^ current[7];
                next[3] = current[3] ^ current[4] ^ current[6] ^ current[7];
                next[4] = current[3] ^ current[4] ^ current[5];
                next[5] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[6] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[7] = current[2] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[5];
                next[1] = current[6];
                next[2] = current[7];
                next[3] = current[5] ^ current[8];
                next[4] = current[6] ^ current[9];
                next[5] = current[0] ^ current[7];
                next[6] = current[1] ^ current[8];
                next[7] = current[2] ^ current[9];
                next[8] = current[3];
                next[9] = current[4];
            end
            default: begin // default to GF(2^8)
                next[0] = current[3] ^ current[7];
                next[1] = current[4];
                next[2] = current[3] ^ current[5] ^ current[7];
                next[3] = current[3] ^ current[4] ^ current[6] ^ current[7];
                next[4] = current[3] ^ current[4] ^ current[5];
                next[5] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[6] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[7] = current[2] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^6
// ==============================================
module gf_mul_alpha_6 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[5];
                next[2] = current[1] ^ current[2];
                next[3] = current[2] ^ current[3];
                next[4] = current[3] ^ current[4];
                next[5] = current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[2] ^ current[6] ^ current[7];
                next[1] = current[3] ^ current[7];
                next[2] = current[2] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[5] ^ current[6];
                next[4] = current[2] ^ current[3] ^ current[4];
                next[5] = current[3] ^ current[4] ^ current[5];
                next[6] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[7] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[4];
                next[1] = current[5];
                next[2] = current[6];
                next[3] = current[4] ^ current[7];
                next[4] = current[5] ^ current[8];
                next[5] = current[6] ^ current[9];
                next[6] = current[0] ^ current[7];
                next[7] = current[1] ^ current[8];
                next[8] = current[2] ^ current[9];
                next[9] = current[3];
            end
            default: begin // default to GF(2^8)
                next[0] = current[2] ^ current[6] ^ current[7];
                next[1] = current[3] ^ current[7];
                next[2] = current[2] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[5] ^ current[6];
                next[4] = current[2] ^ current[3] ^ current[4];
                next[5] = current[3] ^ current[4] ^ current[5];
                next[6] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[7] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^7
// ==============================================
module gf_mul_alpha_7 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[4] ^ current[5];
                next[1] = current[0] ^ current[4];
                next[2] = current[0] ^ current[1] ^ current[5];
                next[3] = current[1] ^ current[2];
                next[4] = current[2] ^ current[3];
                next[5] = current[3] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[1] = current[2] ^ current[6] ^ current[7];
                next[2] = current[1] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[4] ^ current[5];
                next[4] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[5] = current[2] ^ current[3] ^ current[4];
                next[6] = current[3] ^ current[4] ^ current[5];
                next[7] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[3];
                next[1] = current[4];
                next[2] = current[5];
                next[3] = current[3] ^ current[6];
                next[4] = current[4] ^ current[7];
                next[5] = current[5] ^ current[8];
                next[6] = current[6] ^ current[9];
                next[7] = current[0] ^ current[7];
                next[8] = current[1] ^ current[8];
                next[9] = current[2] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[1] = current[2] ^ current[6] ^ current[7];
                next[2] = current[1] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[4] ^ current[5];
                next[4] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[5] = current[2] ^ current[3] ^ current[4];
                next[6] = current[3] ^ current[4] ^ current[5];
                next[7] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^8
// ==============================================
module gf_mul_alpha_8 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[3] ^ current[4];
                next[1] = current[3] ^ current[5];
                next[2] = current[0] ^ current[4];
                next[3] = current[0] ^ current[1] ^ current[5];
                next[4] = current[1] ^ current[2];
                next[5] = current[2] ^ current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[1] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[2] = current[0] ^ current[2] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[3] ^ current[4];
                next[4] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[5] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[6] = current[2] ^ current[3] ^ current[4];
                next[7] = current[3] ^ current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[9];
                next[1] = current[3];
                next[2] = current[4];
                next[3] = current[2] ^ current[5] ^ current[9];
                next[4] = current[3] ^ current[6];
                next[5] = current[4] ^ current[7];
                next[6] = current[5] ^ current[8];
                next[7] = current[6] ^ current[9];
                next[8] = current[0] ^ current[7];
                next[9] = current[1] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[1] = current[1] ^ current[5] ^ current[6] ^ current[7];
                next[2] = current[0] ^ current[2] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[3] ^ current[4];
                next[4] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[5] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[6] = current[2] ^ current[3] ^ current[4];
                next[7] = current[3] ^ current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^9
// ==============================================
module gf_mul_alpha_9 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[2] ^ current[3];
                next[1] = current[2] ^ current[4];
                next[2] = current[3] ^ current[5];
                next[3] = current[0] ^ current[4];
                next[4] = current[0] ^ current[1] ^ current[5];
                next[5] = current[1] ^ current[2];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[3] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[2] = current[1] ^ current[3] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[2] ^ current[3] ^ current[7];
                next[4] = current[0] ^ current[1] ^ current[5];
                next[5] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[6] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[7] = current[2] ^ current[3] ^ current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[1] ^ current[8];
                next[1] = current[2] ^ current[9];
                next[2] = current[3];
                next[3] = current[1] ^ current[4] ^ current[8];
                next[4] = current[2] ^ current[5] ^ current[9];
                next[5] = current[3] ^ current[6];
                next[6] = current[4] ^ current[7];
                next[7] = current[5] ^ current[8];
                next[8] = current[6] ^ current[9];
                next[9] = current[0] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[3] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[4] ^ current[5] ^ current[6];
                next[2] = current[1] ^ current[3] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[2] ^ current[3] ^ current[7];
                next[4] = current[0] ^ current[1] ^ current[5];
                next[5] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[6] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[7] = current[2] ^ current[3] ^ current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^10
// ==============================================
module gf_mul_alpha_10 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[2];
                next[1] = current[1] ^ current[3];
                next[2] = current[2] ^ current[4];
                next[3] = current[3] ^ current[5];
                next[4] = current[0] ^ current[4];
                next[5] = current[0] ^ current[1] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[2] ^ current[3] ^ current[4];
                next[1] = current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[4] ^ current[7];
                next[5] = current[0] ^ current[1] ^ current[5];
                next[6] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[7] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[7];
                next[1] = current[1] ^ current[8];
                next[2] = current[2] ^ current[9];
                next[3] = current[0] ^ current[3] ^ current[7];
                next[4] = current[1] ^ current[4] ^ current[8];
                next[5] = current[2] ^ current[5] ^ current[9];
                next[6] = current[3] ^ current[6];
                next[7] = current[4] ^ current[7];
                next[8] = current[5] ^ current[8];
                next[9] = current[6] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[2] ^ current[3] ^ current[4];
                next[1] = current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[4] ^ current[7];
                next[5] = current[0] ^ current[1] ^ current[5];
                next[6] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[7] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^11
// ==============================================
module gf_mul_alpha_11 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[5];
                next[1] = current[0] ^ current[2] ^ current[5];
                next[2] = current[1] ^ current[3];
                next[3] = current[2] ^ current[4];
                next[4] = current[3] ^ current[5];
                next[5] = current[0] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[1] = current[2] ^ current[3] ^ current[4];
                next[2] = current[1] ^ current[2] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[3] ^ current[6];
                next[5] = current[0] ^ current[4] ^ current[7];
                next[6] = current[0] ^ current[1] ^ current[5];
                next[7] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[6] ^ current[9];
                next[1] = current[0] ^ current[7];
                next[2] = current[1] ^ current[8];
                next[3] = current[2] ^ current[6];
                next[4] = current[0] ^ current[3] ^ current[7];
                next[5] = current[1] ^ current[4] ^ current[8];
                next[6] = current[2] ^ current[5] ^ current[9];
                next[7] = current[3] ^ current[6];
                next[8] = current[4] ^ current[7];
                next[9] = current[5] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[1] = current[2] ^ current[3] ^ current[4];
                next[2] = current[1] ^ current[2] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[3] ^ current[6];
                next[5] = current[0] ^ current[4] ^ current[7];
                next[6] = current[0] ^ current[1] ^ current[5];
                next[7] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^12
// ==============================================
module gf_mul_alpha_12 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[4];
                next[1] = current[1] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[2] ^ current[5];
                next[3] = current[1] ^ current[3];
                next[4] = current[2] ^ current[4];
                next[5] = current[3] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[1] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[3] ^ current[4] ^ current[6];
                next[3] = current[0] ^ current[4] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[2] ^ current[5] ^ current[7];
                next[5] = current[3] ^ current[6];
                next[6] = current[0] ^ current[4] ^ current[7];
                next[7] = current[0] ^ current[1] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[5] ^ current[8];
                next[1] = current[6] ^ current[9];
                next[2] = current[0] ^ current[7];
                next[3] = current[1] ^ current[5];
                next[4] = current[2] ^ current[6];
                next[5] = current[0] ^ current[3] ^ current[7];
                next[6] = current[1] ^ current[4] ^ current[8];
                next[7] = current[2] ^ current[5] ^ current[9];
                next[8] = current[3] ^ current[6];
                next[9] = current[4] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[1] = current[1] ^ current[2] ^ current[3] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[3] ^ current[4] ^ current[6];
                next[3] = current[0] ^ current[4] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[2] ^ current[5] ^ current[7];
                next[5] = current[3] ^ current[6];
                next[6] = current[0] ^ current[4] ^ current[7];
                next[7] = current[0] ^ current[1] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^13
// ==============================================
module gf_mul_alpha_13 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[3] ^ current[5];
                next[1] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[2] = current[1] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[2] ^ current[5];
                next[4] = current[1] ^ current[3];
                next[5] = current[2] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[5] ^ current[7];
                next[3] = current[3] ^ current[4] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[5] = current[2] ^ current[5] ^ current[7];
                next[6] = current[3] ^ current[6];
                next[7] = current[0] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[4] ^ current[7];
                next[1] = current[5] ^ current[8];
                next[2] = current[6] ^ current[9];
                next[3] = current[0] ^ current[4];
                next[4] = current[1] ^ current[5];
                next[5] = current[2] ^ current[6];
                next[6] = current[0] ^ current[3] ^ current[7];
                next[7] = current[1] ^ current[4] ^ current[8];
                next[8] = current[2] ^ current[5] ^ current[9];
                next[9] = current[3] ^ current[6];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[2] ^ current[6];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[5] ^ current[7];
                next[3] = current[3] ^ current[4] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[5] = current[2] ^ current[5] ^ current[7];
                next[6] = current[3] ^ current[6];
                next[7] = current[0] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^14
// ==============================================
module gf_mul_alpha_14 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[2] ^ current[4];
                next[1] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[3] = current[1] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[2] ^ current[5];
                next[5] = current[1] ^ current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[4] ^ current[7];
                next[1] = current[0] ^ current[1] ^ current[5];
                next[2] = current[1] ^ current[2] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[5] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[6] = current[2] ^ current[5] ^ current[7];
                next[7] = current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[3] ^ current[6];
                next[1] = current[4] ^ current[7];
                next[2] = current[5] ^ current[8];
                next[3] = current[3] ^ current[9];
                next[4] = current[0] ^ current[4];
                next[5] = current[1] ^ current[5];
                next[6] = current[2] ^ current[6];
                next[7] = current[0] ^ current[3] ^ current[7];
                next[8] = current[1] ^ current[4] ^ current[8];
                next[9] = current[2] ^ current[5] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[4] ^ current[7];
                next[1] = current[0] ^ current[1] ^ current[5];
                next[2] = current[1] ^ current[2] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[5] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[6] = current[2] ^ current[5] ^ current[7];
                next[7] = current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^15
// ==============================================
module gf_mul_alpha_15 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[3];
                next[1] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[2] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[1] ^ current[4] ^ current[5];
                next[5] = current[0] ^ current[2] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[3] ^ current[6];
                next[1] = current[0] ^ current[4] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[7];
                next[4] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[5] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[6] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[7] = current[2] ^ current[5] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[5] ^ current[9];
                next[1] = current[3] ^ current[6];
                next[2] = current[4] ^ current[7];
                next[3] = current[2] ^ current[8] ^ current[9];
                next[4] = current[3] ^ current[9];
                next[5] = current[0] ^ current[4];
                next[6] = current[1] ^ current[5];
                next[7] = current[2] ^ current[6];
                next[8] = current[0] ^ current[3] ^ current[7];
                next[9] = current[1] ^ current[4] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[3] ^ current[6];
                next[1] = current[0] ^ current[4] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[7];
                next[4] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[5] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[6] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[7] = current[2] ^ current[5] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^16
// ==============================================
module gf_mul_alpha_16 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[2] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[2] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[3] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[5] = current[1] ^ current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[2] ^ current[5] ^ current[7];
                next[1] = current[3] ^ current[6];
                next[2] = current[0] ^ current[2] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[6] ^ current[7];
                next[4] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[5] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[6] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[7] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[1] ^ current[4] ^ current[8];
                next[1] = current[2] ^ current[5] ^ current[9];
                next[2] = current[3] ^ current[6];
                next[3] = current[1] ^ current[7] ^ current[8];
                next[4] = current[2] ^ current[8] ^ current[9];
                next[5] = current[3] ^ current[9];
                next[6] = current[0] ^ current[4];
                next[7] = current[1] ^ current[5];
                next[8] = current[2] ^ current[6];
                next[9] = current[0] ^ current[3] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[2] ^ current[5] ^ current[7];
                next[1] = current[3] ^ current[6];
                next[2] = current[0] ^ current[2] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[6] ^ current[7];
                next[4] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[5] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[6] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[7] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^17
// ==============================================
module gf_mul_alpha_17 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[2] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[3] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[4] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[5] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[1] = current[2] ^ current[5] ^ current[7];
                next[2] = current[1] ^ current[3] ^ current[4] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[5] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[6] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[7] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[3] ^ current[7];
                next[1] = current[1] ^ current[4] ^ current[8];
                next[2] = current[2] ^ current[5] ^ current[9];
                next[3] = current[0] ^ current[6] ^ current[7];
                next[4] = current[1] ^ current[7] ^ current[8];
                next[5] = current[2] ^ current[8] ^ current[9];
                next[6] = current[3] ^ current[9];
                next[7] = current[0] ^ current[4];
                next[8] = current[1] ^ current[5];
                next[9] = current[2] ^ current[6];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[1] = current[2] ^ current[5] ^ current[7];
                next[2] = current[1] ^ current[3] ^ current[4] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[5] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[5] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[6] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[7] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^18
// ==============================================
module gf_mul_alpha_18 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[3];
                next[2] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[4] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[5] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[1] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[6];
                next[3] = current[0] ^ current[1] ^ current[4] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[2] ^ current[3];
                next[5] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[6] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[7] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[6];
                next[1] = current[0] ^ current[3] ^ current[7];
                next[2] = current[1] ^ current[4] ^ current[8];
                next[3] = current[5] ^ current[6] ^ current[9];
                next[4] = current[0] ^ current[6] ^ current[7];
                next[5] = current[1] ^ current[7] ^ current[8];
                next[6] = current[2] ^ current[8] ^ current[9];
                next[7] = current[3] ^ current[9];
                next[8] = current[0] ^ current[4];
                next[9] = current[1] ^ current[5];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[1] = current[1] ^ current[4] ^ current[6] ^ current[7];
                next[2] = current[0] ^ current[2] ^ current[3] ^ current[6];
                next[3] = current[0] ^ current[1] ^ current[4] ^ current[5] ^ current[6];
                next[4] = current[1] ^ current[2] ^ current[3];
                next[5] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[6] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[7] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^19
// ==============================================
module gf_mul_alpha_19 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[2] ^ current[3] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[2];
                next[2] = current[0] ^ current[1] ^ current[3];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[4] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[5] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[1] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[2] = current[1] ^ current[2] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[1] ^ current[2];
                next[5] = current[1] ^ current[2] ^ current[3];
                next[6] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[7] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[1] ^ current[5];
                next[1] = current[2] ^ current[6];
                next[2] = current[0] ^ current[3] ^ current[7];
                next[3] = current[4] ^ current[5] ^ current[8];
                next[4] = current[5] ^ current[6] ^ current[9];
                next[5] = current[0] ^ current[6] ^ current[7];
                next[6] = current[1] ^ current[7] ^ current[8];
                next[7] = current[2] ^ current[8] ^ current[9];
                next[8] = current[3] ^ current[9];
                next[9] = current[0] ^ current[4];
            end
            default: begin // default to GF(2^8)
                next[0] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[1] = current[0] ^ current[3] ^ current[5] ^ current[6] ^ current[7];
                next[2] = current[1] ^ current[2] ^ current[5] ^ current[7];
                next[3] = current[0] ^ current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[1] ^ current[2];
                next[5] = current[1] ^ current[2] ^ current[3];
                next[6] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[7] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^20
// ==============================================
module gf_mul_alpha_20 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[2] ^ current[3] ^ current[4];
                next[1] = current[1] ^ current[5];
                next[2] = current[0] ^ current[2];
                next[3] = current[0] ^ current[1] ^ current[3];
                next[4] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[5] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[1] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[2] = current[0] ^ current[1] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[4] ^ current[7];
                next[4] = current[0] ^ current[1];
                next[5] = current[0] ^ current[1] ^ current[2];
                next[6] = current[1] ^ current[2] ^ current[3];
                next[7] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[4];
                next[1] = current[1] ^ current[5];
                next[2] = current[2] ^ current[6];
                next[3] = current[3] ^ current[4] ^ current[7];
                next[4] = current[4] ^ current[5] ^ current[8];
                next[5] = current[5] ^ current[6] ^ current[9];
                next[6] = current[0] ^ current[6] ^ current[7];
                next[7] = current[1] ^ current[7] ^ current[8];
                next[8] = current[2] ^ current[8] ^ current[9];
                next[9] = current[3] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[1] = current[2] ^ current[4] ^ current[5] ^ current[6];
                next[2] = current[0] ^ current[1] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[4] ^ current[7];
                next[4] = current[0] ^ current[1];
                next[5] = current[0] ^ current[1] ^ current[2];
                next[6] = current[1] ^ current[2] ^ current[3];
                next[7] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^21
// ==============================================
module gf_mul_alpha_21 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[5];
                next[1] = current[0] ^ current[4] ^ current[5];
                next[2] = current[1] ^ current[5];
                next[3] = current[0] ^ current[2];
                next[4] = current[0] ^ current[1] ^ current[3];
                next[5] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[1] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[3] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[7];
                next[5] = current[0] ^ current[1];
                next[6] = current[0] ^ current[1] ^ current[2];
                next[7] = current[1] ^ current[2] ^ current[3];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[3] ^ current[9];
                next[1] = current[0] ^ current[4];
                next[2] = current[1] ^ current[5];
                next[3] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[4] = current[3] ^ current[4] ^ current[7];
                next[5] = current[4] ^ current[5] ^ current[8];
                next[6] = current[5] ^ current[6] ^ current[9];
                next[7] = current[0] ^ current[6] ^ current[7];
                next[8] = current[1] ^ current[7] ^ current[8];
                next[9] = current[2] ^ current[8] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[1] = current[1] ^ current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[3] ^ current[5] ^ current[6];
                next[3] = current[1] ^ current[2] ^ current[3] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[7];
                next[5] = current[0] ^ current[1];
                next[6] = current[0] ^ current[1] ^ current[2];
                next[7] = current[1] ^ current[2] ^ current[3];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^22
// ==============================================
module gf_mul_alpha_22 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[2] ^ current[4];
                next[1] = current[3] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[4] ^ current[5];
                next[3] = current[1] ^ current[5];
                next[4] = current[0] ^ current[2];
                next[5] = current[0] ^ current[1] ^ current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[3];
                next[1] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[2] = current[2] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[5] ^ current[6];
                next[4] = current[6] ^ current[7];
                next[5] = current[0] ^ current[7];
                next[6] = current[0] ^ current[1];
                next[7] = current[0] ^ current[1] ^ current[2];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[8] ^ current[9];
                next[1] = current[3] ^ current[9];
                next[2] = current[0] ^ current[4];
                next[3] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[4] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[5] = current[3] ^ current[4] ^ current[7];
                next[6] = current[4] ^ current[5] ^ current[8];
                next[7] = current[5] ^ current[6] ^ current[9];
                next[8] = current[0] ^ current[6] ^ current[7];
                next[9] = current[1] ^ current[7] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[3];
                next[1] = current[0] ^ current[2] ^ current[3] ^ current[4];
                next[2] = current[2] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[5] ^ current[6];
                next[4] = current[6] ^ current[7];
                next[5] = current[0] ^ current[7];
                next[6] = current[0] ^ current[1];
                next[7] = current[0] ^ current[1] ^ current[2];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^23
// ==============================================
module gf_mul_alpha_23 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[3];
                next[1] = current[2] ^ current[3] ^ current[4];
                next[2] = current[3] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[4] ^ current[5];
                next[4] = current[1] ^ current[5];
                next[5] = current[0] ^ current[2];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[2];
                next[1] = current[1] ^ current[2] ^ current[3];
                next[2] = current[1] ^ current[3] ^ current[4];
                next[3] = current[0] ^ current[1] ^ current[4] ^ current[5];
                next[4] = current[5] ^ current[6];
                next[5] = current[6] ^ current[7];
                next[6] = current[0] ^ current[7];
                next[7] = current[0] ^ current[1];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[1] ^ current[7] ^ current[8];
                next[1] = current[2] ^ current[8] ^ current[9];
                next[2] = current[3] ^ current[9];
                next[3] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[4] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[5] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[6] = current[3] ^ current[4] ^ current[7];
                next[7] = current[4] ^ current[5] ^ current[8];
                next[8] = current[5] ^ current[6] ^ current[9];
                next[9] = current[0] ^ current[6] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[2];
                next[1] = current[1] ^ current[2] ^ current[3];
                next[2] = current[1] ^ current[3] ^ current[4];
                next[3] = current[0] ^ current[1] ^ current[4] ^ current[5];
                next[4] = current[5] ^ current[6];
                next[5] = current[6] ^ current[7];
                next[6] = current[0] ^ current[7];
                next[7] = current[0] ^ current[1];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^24
// ==============================================
module gf_mul_alpha_24 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[2];
                next[1] = current[1] ^ current[2] ^ current[3];
                next[2] = current[2] ^ current[3] ^ current[4];
                next[3] = current[3] ^ current[4] ^ current[5];
                next[4] = current[0] ^ current[4] ^ current[5];
                next[5] = current[1] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[1];
                next[1] = current[0] ^ current[1] ^ current[2];
                next[2] = current[0] ^ current[2] ^ current[3];
                next[3] = current[0] ^ current[3] ^ current[4];
                next[4] = current[4] ^ current[5];
                next[5] = current[5] ^ current[6];
                next[6] = current[6] ^ current[7];
                next[7] = current[0] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[6] ^ current[7];
                next[1] = current[1] ^ current[7] ^ current[8];
                next[2] = current[2] ^ current[8] ^ current[9];
                next[3] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[4] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[5] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[6] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[7] = current[3] ^ current[4] ^ current[7];
                next[8] = current[4] ^ current[5] ^ current[8];
                next[9] = current[5] ^ current[6] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[1];
                next[1] = current[0] ^ current[1] ^ current[2];
                next[2] = current[0] ^ current[2] ^ current[3];
                next[3] = current[0] ^ current[3] ^ current[4];
                next[4] = current[4] ^ current[5];
                next[5] = current[5] ^ current[6];
                next[6] = current[6] ^ current[7];
                next[7] = current[0] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^25
// ==============================================
module gf_mul_alpha_25 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[2] = current[1] ^ current[2] ^ current[3];
                next[3] = current[2] ^ current[3] ^ current[4];
                next[4] = current[3] ^ current[4] ^ current[5];
                next[5] = current[0] ^ current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[7];
                next[1] = current[0] ^ current[1];
                next[2] = current[1] ^ current[2] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[7];
                next[4] = current[3] ^ current[4] ^ current[7];
                next[5] = current[4] ^ current[5];
                next[6] = current[5] ^ current[6];
                next[7] = current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[5] ^ current[6] ^ current[9];
                next[1] = current[0] ^ current[6] ^ current[7];
                next[2] = current[1] ^ current[7] ^ current[8];
                next[3] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[4] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[5] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[6] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[7] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[8] = current[3] ^ current[4] ^ current[7];
                next[9] = current[4] ^ current[5] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[7];
                next[1] = current[0] ^ current[1];
                next[2] = current[1] ^ current[2] ^ current[7];
                next[3] = current[2] ^ current[3] ^ current[7];
                next[4] = current[3] ^ current[4] ^ current[7];
                next[5] = current[4] ^ current[5];
                next[6] = current[5] ^ current[6];
                next[7] = current[6] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^26
// ==============================================
module gf_mul_alpha_26 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[1] ^ current[4];
                next[2] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[3] = current[1] ^ current[2] ^ current[3];
                next[4] = current[2] ^ current[3] ^ current[4];
                next[5] = current[3] ^ current[4] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[6] ^ current[7];
                next[1] = current[0] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[6] ^ current[7];
                next[3] = current[1] ^ current[2] ^ current[6];
                next[4] = current[2] ^ current[3] ^ current[6];
                next[5] = current[3] ^ current[4] ^ current[7];
                next[6] = current[4] ^ current[5];
                next[7] = current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[4] ^ current[5] ^ current[8];
                next[1] = current[5] ^ current[6] ^ current[9];
                next[2] = current[0] ^ current[6] ^ current[7];
                next[3] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[4] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[5] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[6] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[7] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[8] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[9] = current[3] ^ current[4] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[6] ^ current[7];
                next[1] = current[0] ^ current[7];
                next[2] = current[0] ^ current[1] ^ current[6] ^ current[7];
                next[3] = current[1] ^ current[2] ^ current[6];
                next[4] = current[2] ^ current[3] ^ current[6];
                next[5] = current[3] ^ current[4] ^ current[7];
                next[6] = current[4] ^ current[5];
                next[7] = current[5] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^27
// ==============================================
module gf_mul_alpha_27 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[3] ^ current[4] ^ current[5];
                next[1] = current[0] ^ current[3];
                next[2] = current[0] ^ current[1] ^ current[4];
                next[3] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[4] = current[1] ^ current[2] ^ current[3];
                next[5] = current[2] ^ current[3] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[5] ^ current[6];
                next[1] = current[6] ^ current[7];
                next[2] = current[0] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[5] ^ current[7];
                next[4] = current[1] ^ current[2] ^ current[5];
                next[5] = current[2] ^ current[3] ^ current[6];
                next[6] = current[3] ^ current[4] ^ current[7];
                next[7] = current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[3] ^ current[4] ^ current[7];
                next[1] = current[4] ^ current[5] ^ current[8];
                next[2] = current[5] ^ current[6] ^ current[9];
                next[3] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[4] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[5] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[6] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[7] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[8] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[9] = current[2] ^ current[3] ^ current[6] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[5] ^ current[6];
                next[1] = current[6] ^ current[7];
                next[2] = current[0] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[1] ^ current[5] ^ current[7];
                next[4] = current[1] ^ current[2] ^ current[5];
                next[5] = current[2] ^ current[3] ^ current[6];
                next[6] = current[3] ^ current[4] ^ current[7];
                next[7] = current[4] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^28
// ==============================================
module gf_mul_alpha_28 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[2] ^ current[3] ^ current[4];
                next[1] = current[2] ^ current[5];
                next[2] = current[0] ^ current[3];
                next[3] = current[0] ^ current[1] ^ current[4];
                next[4] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[5] = current[1] ^ current[2] ^ current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[4] ^ current[5];
                next[1] = current[5] ^ current[6];
                next[2] = current[4] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[4] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[5] = current[1] ^ current[2] ^ current[5];
                next[6] = current[2] ^ current[3] ^ current[6];
                next[7] = current[3] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[1] = current[3] ^ current[4] ^ current[7];
                next[2] = current[4] ^ current[5] ^ current[8];
                next[3] = current[2] ^ current[3] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[5] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[6] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[7] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[8] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[9] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[4] ^ current[5];
                next[1] = current[5] ^ current[6];
                next[2] = current[4] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[4] ^ current[6] ^ current[7];
                next[4] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[5] = current[1] ^ current[2] ^ current[5];
                next[6] = current[2] ^ current[3] ^ current[6];
                next[7] = current[3] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^29
// ==============================================
module gf_mul_alpha_29 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[1] ^ current[2] ^ current[3];
                next[1] = current[1] ^ current[4];
                next[2] = current[2] ^ current[5];
                next[3] = current[0] ^ current[3];
                next[4] = current[0] ^ current[1] ^ current[4];
                next[5] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[3] ^ current[4] ^ current[7];
                next[1] = current[4] ^ current[5];
                next[2] = current[3] ^ current[4] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[3] ^ current[5] ^ current[6];
                next[4] = current[0] ^ current[3] ^ current[6];
                next[5] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[6] = current[1] ^ current[2] ^ current[5];
                next[7] = current[2] ^ current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[1] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[2] = current[3] ^ current[4] ^ current[7];
                next[3] = current[1] ^ current[2] ^ current[4] ^ current[9];
                next[4] = current[2] ^ current[3] ^ current[5];
                next[5] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[6] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[7] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[8] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[9] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[3] ^ current[4] ^ current[7];
                next[1] = current[4] ^ current[5];
                next[2] = current[3] ^ current[4] ^ current[5] ^ current[6] ^ current[7];
                next[3] = current[3] ^ current[5] ^ current[6];
                next[4] = current[0] ^ current[3] ^ current[6];
                next[5] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[6] = current[1] ^ current[2] ^ current[5];
                next[7] = current[2] ^ current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^30
// ==============================================
module gf_mul_alpha_30 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[2] ^ current[5];
                next[1] = current[0] ^ current[3] ^ current[5];
                next[2] = current[1] ^ current[4];
                next[3] = current[2] ^ current[5];
                next[4] = current[0] ^ current[3];
                next[5] = current[0] ^ current[1] ^ current[4];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[2] ^ current[3] ^ current[6];
                next[1] = current[3] ^ current[4] ^ current[7];
                next[2] = current[2] ^ current[3] ^ current[4] ^ current[5] ^ current[6];
                next[3] = current[2] ^ current[4] ^ current[5] ^ current[7];
                next[4] = current[2] ^ current[5];
                next[5] = current[0] ^ current[3] ^ current[6];
                next[6] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[7] = current[1] ^ current[2] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[1] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[2] = current[2] ^ current[3] ^ current[6] ^ current[9];
                next[3] = current[0] ^ current[1] ^ current[3] ^ current[8];
                next[4] = current[1] ^ current[2] ^ current[4] ^ current[9];
                next[5] = current[2] ^ current[3] ^ current[5];
                next[6] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[7] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[8] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[9] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
            end
            default: begin // default to GF(2^8)
                next[0] = current[2] ^ current[3] ^ current[6];
                next[1] = current[3] ^ current[4] ^ current[7];
                next[2] = current[2] ^ current[3] ^ current[4] ^ current[5] ^ current[6];
                next[3] = current[2] ^ current[4] ^ current[5] ^ current[7];
                next[4] = current[2] ^ current[5];
                next[5] = current[0] ^ current[3] ^ current[6];
                next[6] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[7] = current[1] ^ current[2] ^ current[5];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^31
// ==============================================
module gf_mul_alpha_31 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[1] ^ current[4];
                next[1] = current[2] ^ current[4] ^ current[5];
                next[2] = current[0] ^ current[3] ^ current[5];
                next[3] = current[1] ^ current[4];
                next[4] = current[2] ^ current[5];
                next[5] = current[0] ^ current[3];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[5];
                next[1] = current[2] ^ current[3] ^ current[6];
                next[2] = current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[1] ^ current[3] ^ current[4] ^ current[6];
                next[4] = current[1] ^ current[4] ^ current[7];
                next[5] = current[2] ^ current[5];
                next[6] = current[0] ^ current[3] ^ current[6];
                next[7] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[1] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[2] = current[1] ^ current[2] ^ current[5] ^ current[8] ^ current[9];
                next[3] = current[0] ^ current[2] ^ current[7];
                next[4] = current[0] ^ current[1] ^ current[3] ^ current[8];
                next[5] = current[1] ^ current[2] ^ current[4] ^ current[9];
                next[6] = current[2] ^ current[3] ^ current[5];
                next[7] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[8] = current[1] ^ current[4] ^ current[5] ^ current[7];
                next[9] = current[2] ^ current[5] ^ current[6] ^ current[8];
            end
            default: begin // default to GF(2^8)
                next[0] = current[1] ^ current[2] ^ current[5];
                next[1] = current[2] ^ current[3] ^ current[6];
                next[2] = current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[5] ^ current[7];
                next[3] = current[1] ^ current[3] ^ current[4] ^ current[6];
                next[4] = current[1] ^ current[4] ^ current[7];
                next[5] = current[2] ^ current[5];
                next[6] = current[0] ^ current[3] ^ current[6];
                next[7] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

// ==============================================
// GF multiply by alpha^32
// ==============================================
module gf_mul_alpha_32 (
    input  [9:0] current,
    input  [1:0] m,
    output reg [9:0] next
);
    always @* begin
        case (m)
            2'b01: begin // GF(2^6)
                next[0] = current[0] ^ current[3];
                next[1] = current[1] ^ current[3] ^ current[4];
                next[2] = current[2] ^ current[4] ^ current[5];
                next[3] = current[0] ^ current[3] ^ current[5];
                next[4] = current[1] ^ current[4];
                next[5] = current[2] ^ current[5];
                next[6] = 1'b0;
                next[7] = 1'b0;
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b10: begin // GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[1] = current[1] ^ current[2] ^ current[5];
                next[2] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[2] ^ current[3] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[6] ^ current[7];
                next[5] = current[1] ^ current[4] ^ current[7];
                next[6] = current[2] ^ current[5];
                next[7] = current[0] ^ current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
            2'b11: begin // GF(2^10)
                next[0] = current[2] ^ current[5] ^ current[6] ^ current[8];
                next[1] = current[0] ^ current[3] ^ current[6] ^ current[7] ^ current[9];
                next[2] = current[0] ^ current[1] ^ current[4] ^ current[7] ^ current[8];
                next[3] = current[1] ^ current[6] ^ current[9];
                next[4] = current[0] ^ current[2] ^ current[7];
                next[5] = current[0] ^ current[1] ^ current[3] ^ current[8];
                next[6] = current[1] ^ current[2] ^ current[4] ^ current[9];
                next[7] = current[2] ^ current[3] ^ current[5];
                next[8] = current[0] ^ current[3] ^ current[4] ^ current[6];
                next[9] = current[1] ^ current[4] ^ current[5] ^ current[7];
            end
            default: begin // default to GF(2^8)
                next[0] = current[0] ^ current[1] ^ current[4] ^ current[7];
                next[1] = current[1] ^ current[2] ^ current[5];
                next[2] = current[0] ^ current[1] ^ current[2] ^ current[3] ^ current[4] ^ current[6] ^ current[7];
                next[3] = current[0] ^ current[2] ^ current[3] ^ current[5];
                next[4] = current[0] ^ current[3] ^ current[6] ^ current[7];
                next[5] = current[1] ^ current[4] ^ current[7];
                next[6] = current[2] ^ current[5];
                next[7] = current[0] ^ current[3] ^ current[6];
                next[8] = 1'b0;
                next[9] = 1'b0;
            end
        endcase
    end
endmodule

