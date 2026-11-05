`timescale 1ns/1ps
module decoder_vectors;
  bit clk = 1, rstn, mode, set;
  logic [1:0] code;
  logic [63:0] idata;
  wire ready, finish;
  wire [9:0] odata;
  bch dut(.*);
  always #5 clk = !clk;

  logic [63:0] inputs[$];
  logic [9:0] answers[$];
  int sent, received;
  always @(negedge clk) begin
    if (ready && sent < inputs.size()) begin
      idata = inputs[sent];
      sent++;
    end
    if (rstn && finish) begin
      if (received >= answers.size() || odata !== answers[received])
        $fatal(1, "Result %0d: got %0d, expected %0d", received, odata, answers[received]);
      received++;
    end
  end

  initial begin : run
    int pattern, fd, words;
    logic [63:0] word;
    logic [9:0] position;
    if (!$value$plusargs("PATTERN=%d", pattern)) pattern = 100;
    if (!(pattern inside {100,200,300,400,500,600})) $fatal(1,"Unknown fixture");
    code = 2'(((pattern / 100 - 1) % 3) + 1);
    mode = pattern >= 400;
    fd = $fopen($sformatf("sim/vectors/p%0d.txt",pattern),"r");
    if (!fd) $fatal(1,"Input file unavailable");
    while ($fscanf(fd,"%b",word)==1) inputs.push_back(word);
    $fclose(fd);
    fd = $fopen($sformatf("sim/vectors/p%0da.txt",pattern),"r");
    if (!fd) $fatal(1,"Answer file unavailable");
    while ($fscanf(fd,"%b",position)==1) answers.push_back(position);
    $fclose(fd);
    words = code == 1 ? 8 : code == 2 ? 32 : 128;
    if (!inputs.size() || inputs.size() % words || !answers.size()) $fatal(1,"Incomplete fixture");
    repeat (6) @(negedge clk);
    rstn = 1;
    repeat (5) @(negedge clk);
    repeat (inputs.size()/words) begin
      set = 1;
      @(negedge clk); set = 0;
      wait (finish);
      wait (!finish);
      repeat (4) @(negedge clk);
    end
    if (received != answers.size() || sent != inputs.size())
      $fatal(1,"Consumed %0d/%0d words and %0d/%0d answers",sent,inputs.size(),received,answers.size());
    $display("PASS p%0d: %0d words, %0d error locations",pattern,sent,received);
    $finish;
  end
  initial #1000000 $fatal(1,"Timeout");
endmodule
