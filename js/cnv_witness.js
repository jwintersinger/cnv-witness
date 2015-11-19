function CnvPlotter() {
  var horiz_padding = 50;
  var vert_padding = 30;
  this._M = [vert_padding, horiz_padding, vert_padding, horiz_padding],
      this._W = 1200 - this._M[1] - this._M[3],
      this._H = 500 - this._M[0] - this._M[2];

  var svg = d3.select('#container').html('')
      .append('svg:svg')
      .attr('width', this._W + this._M[1] + this._M[3])
      .attr('height', this._H + this._M[0] + this._M[2]);
  this._container = svg.append('svg:g')  // This container is for padding
                       .attr('transform', 'translate(' + this._M[3] + ',' + this._M[0] + ')')
                       .append('svg:g'); // This container is for zooming

  var self = this;
  svg.call(d3.behavior.zoom().on('zoom', function() {
    self._container.attr('transform', 'translate(' + d3.event.translate + ') scale(' + d3.event.scale + ')');
  }).scaleExtent([1, 100]));

  var max_cn = 5;
  this._compute_chrom_lens();

  this._xscale = d3.scale.linear()
                         .domain([1, this._total_len])
                         .range([0, this._W]);
  this._yscale = d3.scale.linear()
                         .domain([0, max_cn])
                         .range([this._H, 0]);

  var yticks = [];
  for(var i = 0; i <= max_cn; i++) { yticks.push(i); }
  this._yaxis = d3.svg.axis()
                      .scale(this._yscale)
                      .tickValues(yticks)
                      .tickFormat(d3.format('d'))
                      .orient('left');
  this._container.append('svg:g')
           .attr('class', 'axis')
           .call(this._yaxis);
  this._draw_chrom_markers();
}

CnvPlotter.prototype._compute_chrom_lens = function() {
  var chr_lens = [
    ['1', 248956422],
    ['2', 242193529],
    ['3', 198295559],
    ['4', 190214555],
    ['5', 181538259],
    ['6', 170805979],
    ['7', 159345973],
    ['8', 145138636],
    ['9', 138394717],
    ['10', 133797422],
    ['11', 135086622],
    ['12', 133275309],
    ['13', 114364328],
    ['14', 107043718],
    ['15', 101991189],
    ['16', 90338345],
    ['17', 83257441],
    ['18', 80373285],
    ['19', 58617616],
    ['20', 64444167],
    ['21', 46709983],
    ['22', 50818468],
    ['X', 156040895],
    ['Y', 57227415]
  ];

  var self = this;

  this._chr_lens = {};
  chr_lens.forEach(function(pair) {
    var chrom = pair[0], len = pair[1];
    self._chr_lens[chrom] = len;
  });

  var sum = 0;
  self._cum_chr_lens = {};
  for(var i = 0; i < chr_lens.length; i++) {
    var chr = chr_lens[i][0], size = chr_lens[i][1];
    self._cum_chr_lens[chr] = sum;
    sum += size;
  }

  this._total_len = sum;
}


CnvPlotter.prototype._draw_chrom_markers = function() {
  var chroms = this._container.append('svg:g')
                 .attr('class', 'axis')
                 .selectAll('.chrom')
                 .data(Object.keys(this._cum_chr_lens)).enter()
                 .append('svg:g')
                 .attr('class', 'chrom');

  var self = this;
  var get_xpos_line = function(d, i) {
    return self._xscale(self._cum_chr_lens[d] + self._chr_lens[d]);
  };
  var get_xpos_label = function(d, i) {
    return self._xscale(self._cum_chr_lens[d] + 0.5 * self._chr_lens[d]);
  };

  chroms.append('line')
        .attr('x1', get_xpos_line)
        .attr('x2', get_xpos_line)
        .attr('y1', 0)
        .attr('y2', this._H);
  chroms.append('text')
        .attr('text-anchor', 'middle')
        .attr('x', get_xpos_label)
        .attr('y', this._H + this._M[0])
        .text(function(d, i) { return d; });
}

CnvPlotter.prototype._compute_cum_chr_locus = function(chrom, pos) {
  return this._cum_chr_lens[chrom] + pos;
}

CnvPlotter.prototype._pick_colours = function() {
  return {
    vanloo_wedge: '#ff0000',
    mustonen095: '#00ff00',
    peifer: '#0000ff'
  };
}

CnvPlotter.prototype.plot = function(intervals) {
  var rect_height = 20;
  var colours = this._pick_colours();

  var self = this;
  Object.keys(intervals).forEach(function(chrom) {
    Object.keys(intervals[chrom]).forEach(function(cnstate) {
      var comps = cnstate.split(',');
      var major = parseInt(comps[0], 10), minor = parseInt(comps[1], 10);
      intervals[chrom][cnstate].forEach(function(interval) {
        var start = interval[0], end = interval[1], methods = interval[2];

        [major + minor, major].forEach(function(cn) {
          var xstart = self._xscale(self._compute_cum_chr_locus(chrom, start));
          var xend = self._xscale(self._compute_cum_chr_locus(chrom, end));
          var ystart = self._yscale(cn);
          var yoffset = (rect_height * methods.length) / 2;

          self._container.append('svg:g')
                         .attr('transform', 'translate(0,' + (ystart - yoffset) + ')')
                         .selectAll('rect')
                         .data(methods)
                         .enter().append('rect')
                         .attr('x', xstart)
                         .attr('y', function(d, i) { /*console.log([d, chrom, start, end, xstart, xend, ystart, yoffset, ystart - yoffset + i*rect_height]);*/ return i * rect_height; })
                         .attr('fill', function(d, i) { return colours[d] })
                         .attr('width', xend - xstart)
                         .attr('height', rect_height);
        });
      });
    });
  });
}

function main() {
  var plotter = new CnvPlotter();

  d3.json("data/test.json", function(error, json) {
    if(error) return console.warn(error);
    plotter.plot(json);
  });
}

main();
