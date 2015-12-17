function CnvPlotter(iface) {
  this._iface = iface;

  var horiz_padding = 50;
  var vert_padding = 30;
  this._M = [vert_padding, horiz_padding, vert_padding, horiz_padding],
      //this._W = 1200 - this._M[1] - this._M[3],
      this._H = 500 - this._M[0] - this._M[2];

  this._svg = d3.select('#container').html('')
      .append('svg:svg')
      .attr('width', '100%')
      .attr('height', this._H + this._M[0] + this._M[2]);
  this._W = this._svg.node().getBoundingClientRect().width - this._M[1] - this._M[3];
  this._container = this._svg.append('svg:g')  // This container is for padding
                       .attr('transform', 'translate(' + this._M[3] + ',' + this._M[0] + ')')
                       .append('svg:g'); // This container is for zooming

  var self = this;
  this._svg.call(d3.behavior.zoom().on('zoom', function() {
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

CnvPlotter.prototype._create_fills = function(colours, rect_heights) {
  var fills = {
    total: colours,
    minor: {},
  };

  var self = this;
  Object.keys(colours).forEach(function(method) {
    // See http://stackoverflow.com/a/22643745
    var fill_name = method + '_minor';
    self._svg.append('pattern')
       .attr({
         id: fill_name,
         width: 10,
         height: 10,
         patternTransform: 'rotate(45 0 0)',
         patternUnits: 'userSpaceOnUse'
      }).append('line').attr({
        x1: 0,
        y1: 0,
        x2: 0,
        y2: 10,
      }).style({
        'stroke': colours[method],
        'stroke-width': 16
      });
    fills.minor[method] = 'url(#' + fill_name + ')';
  });

  return fills;
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

CnvPlotter.prototype._compute_offsets = function(cn_calls, rect_heights, active_methods) {
  var states = {};

  Object.keys(rect_heights).forEach(function(cntype) {
    var rect_height = rect_heights[cntype];
    Object.keys(cn_calls.intervals[cntype]).forEach(function(chrom) {
      Object.keys(cn_calls.intervals[cntype][chrom]).forEach(function(cnstate) {
        if(!(cnstate in states))
          states[cnstate] = {};
        if(!(cntype in states[cnstate]))
          states[cnstate][cntype] = {};

        cn_calls.intervals[cntype][chrom][cnstate].forEach(function(interval) {
          var start = interval[0], end = interval[1], methods = interval[2];
          methods.forEach(function(method) {
            if(active_methods.hasOwnProperty(method))
              states[cnstate][cntype][method] = true;
          });
        });
      });
    });
  });

  var offsets = {};
  var summed = {};

  Object.keys(states).forEach(function(cnstate) {
    offsets[cnstate] = {};
    // Explicitly list cntype rather than using Object.keys(), as we want to
    // iterate in this exact order.
    var total_cn_offset_sum = 0;
    ['total', 'minor'].forEach(function(cntype) {
      if(!(cntype in states[cnstate]))
        return;

      offsets[cnstate][cntype] = {};
      var methods_for_state = Object.keys(states[cnstate][cntype]).sort();

      for(var i = 0; i < methods_for_state.length; i++) {
        var method = methods_for_state[i];
        offsets[cnstate][cntype][method] = i * rect_heights[cntype];

        if(cntype == 'total')
          total_cn_offset_sum += rect_heights.total;
        else if(cntype == 'minor')
          offsets[cnstate][cntype][method] += total_cn_offset_sum;
        summed[cnstate] = offsets[cnstate][cntype][method] + rect_heights[cntype];
      }
    });
  });

  return {
    indiv: offsets,
    summed: summed
  };
}

CnvPlotter.prototype.plot = function(cn_calls, active_methods) {
  var rect_heights = {
    total: 10,
    minor: 5
  };
  // Object gives O(1) lookup and simpler code.
  var active_meth = {};
  active_methods.forEach(function(method) {
    active_meth[method] = true;
  });
  active_methods = active_meth;

  var offsets = this._compute_offsets(cn_calls, rect_heights, active_methods);
  var colours = this._iface.pick_colours();
  var fills = this._create_fills(colours, rect_heights);

  var self = this;
  Object.keys(rect_heights).forEach(function(cntype) {
    Object.keys(cn_calls.intervals[cntype]).forEach(function(chrom) {
      Object.keys(cn_calls.intervals[cntype][chrom]).forEach(function(cnstate) {
        cn_calls.intervals[cntype][chrom][cnstate].forEach(function(interval) {
          var start = interval[0], end = interval[1], methods = interval[2];

          var meth = [];
          methods.forEach(function(m) {
            if(active_methods.hasOwnProperty(m))
              meth.push(m);
          });
          methods = meth;
          if(methods.length === 0)
            return;

          var xstart = self._xscale(self._compute_cum_chr_locus(chrom, start));
          var xend = self._xscale(self._compute_cum_chr_locus(chrom, end));
          var ystart = self._yscale(cnstate);
          var group_yoffset = offsets.summed[cnstate] / 2;

          self._container.append('svg:g')
                         .attr('transform', 'translate(0,' + (ystart - group_yoffset) + ')')
                         .selectAll('rect')
                         .data(methods)
                         .enter().append('rect')
                         .attr('x', xstart)
                         .attr('y', function(d, i) { return offsets.indiv[cnstate][cntype][d]; })
                         .attr('fill', function(d, i) { return fills[cntype][d] })
                         .attr('width', xend - xstart)
                         .attr('height', rect_heights[cntype]);
        });
      });
    });
  });
}

function Interface(sample_list, metadata) {
  this._index_samples_by_method(sample_list);
  this._fill_sample_selectors(sample_list, metadata);
  this._activate_filters(sample_list);
  this._filter();
}

Interface.prototype.make_methods_filter = function(container, methods, on_change) {
  var colours = this.pick_colours();
  var buttons_container = d3.select(container).html('')
    .append('div')
    .attr('class', 'method-selector btn-group')
    .attr('role', 'group')
    .selectAll('button')
    .data(methods)
    .enter().append('button')
    .attr('class', 'btn btn-primary active method-choice')
    .style('background-color', function(d, i) { return colours[d]; })
    .style('color', 'white')
    .text(function(d, i) { return d; })
    .on('click', function(method) {
      var elem = d3.select(this);
      // Toggle "active" class
      elem.classed('active', !elem.classed('active'));

      var active_methods = d3.select(this.parentNode).selectAll('.method-choice')
                             .filter('.active')
                             .data();
      on_change(active_methods);
    });
}

Interface.prototype._fill_sample_selectors = function(sample_list, metadata) {
  var sampids = Object.keys(sample_list).sort();
  this._update_sample_count(sampids.length);

  var make_human_readable_from_sort_val = function() {
    var prop = parseFloat(d3.select(this).attr('data-sort-value'));
    return prop === -1 ? '&mdash;' : prop.toFixed(3);
  };

  var calc_mean_val = function(sampid, prop) {
    var vals = [];
    if(!(metadata.hasOwnProperty(sampid) && metadata[sampid].hasOwnProperty(prop)))
      return -1;
    Object.keys(metadata[sampid][prop]).forEach(function(method) {
      if(method.lastIndexOf('theta', 0) === 0)
        return;
      vals.push(metadata[sampid][prop][method]);
    });
    return d3.mean(vals);
  };

  // Fill extended selector
  var rows = d3.select('#sample-list-extended tbody').html('')
    .selectAll('tr')
    .data(sampids)
    .enter().append('tr');
  rows.append('td').attr('class', 'sampid').text(function(sampid) { return sampid; });
  rows.append('td').attr('class', 'tumor-type').html(function(sampid) {
    return metadata.hasOwnProperty(sampid) ? metadata[sampid]['tumor_type'] : '&mdash;';
  });
  rows.append('td').attr('class', 'ploidy').attr('data-sort-value', function(sampid) {
    return calc_mean_val(sampid, 'ploidy');
  }).html(make_human_readable_from_sort_val);
  rows.append('td').attr('class', 'purity').attr('data-sort-value', function(sampid) {
    return calc_mean_val(sampid, 'purity');
  }).html(make_human_readable_from_sort_val);
  rows.append('td').attr('class', 'genome-prop').attr('data-sort-value', function(d, i) {
    var proportions = sample_list[d].genome_proportions
    // Return -1 so that datasets without consensus will be sorted below those
    // for which the consensus is actually 0.
    return proportions.hasOwnProperty('consensus') ?  proportions.consensus  : -1;
  }).html(make_human_readable_from_sort_val);
  rows.append('td').attr('class', 'consensus-score').attr('data-sort-value', function(d, i) {
    var proportions = sample_list[d].genome_proportions
    if(!proportions.hasOwnProperty('consensus'))
      return -1;

    // Don't include consensus in mean.
    var cons_prop = proportions.consensus;
    delete proportions.consensus;

    // Don't include THetA values in mean.
    var keys = Object.keys(proportions);
    keys.forEach(function(k) {
      if(k.lastIndexOf('theta', 0) === 0)
        delete proportions[k];
    });

    var mean_prop = d3.mean(Object.keys(proportions).map(function(method) {
      return proportions[method];
    }));
    return cons_prop / mean_prop;
  }).html(make_human_readable_from_sort_val);

  $('#sample-list-extended').stupidtable();

  var self = this;
  var load_sample = function(sampid) {
    var cn_calls_path = sample_list[sampid].cn_calls_path;
    d3.json(cn_calls_path, function(error, cn_calls) {
      var redraw = function() {
        if(error) return console.warn(error);
        var methods = cn_calls.methods;

        var plotter = new CnvPlotter(self);
        plotter.plot(cn_calls, methods);

        self.make_methods_filter('#sample-method-filter', methods, function(active_methods) {
          plotter = new CnvPlotter(self);
          plotter.plot(cn_calls, active_methods);
        });
      };
      redraw();
      d3.select(window).on('resize', redraw);
    });

    d3.select('#sampid').text(sampid);
    if(!metadata.hasOwnProperty(sampid))
      return;
    d3.select('#tumor-type').text(metadata[sampid].tumor_type);

    var sample_meta = {};
    Object.keys(metadata[sampid].ploidy).forEach(function(method) {
      if(!sample_meta.hasOwnProperty(method))
        sample_meta[method] = {};
      sample_meta[method].ploidy = metadata[sampid].ploidy[method];
    });
    Object.keys(metadata[sampid].purity).forEach(function(method) {
      if(!sample_meta.hasOwnProperty(method))
        sample_meta[method] = {};
      sample_meta[method].purity = metadata[sampid].purity[method];
    });

    var rows = d3.select('#sample-metadata tbody').html('')
      .selectAll('tr')
      .data(Object.keys(sample_meta).sort())
      .enter().append('tr');
    rows.append('td').html(function(method) {
      return method;
    });
    rows.append('td').html(function(method) {
      return sample_meta.hasOwnProperty(method) && sample_meta[method].hasOwnProperty('purity') ? sample_meta[method].purity.toFixed(3) : '&mdash;';
    });
    rows.append('td').html(function(method) {
      return sample_meta.hasOwnProperty(method) && sample_meta[method].hasOwnProperty('ploidy') ? sample_meta[method].ploidy.toFixed(3) : '&mdash;';
    });
  };

  d3.selectAll('#sample-list-extended tbody tr').on('click', function(sampid) {
    $('#sample-selector-extended').modal('hide');
    load_sample(sampid);
  });
}

Interface.prototype._update_sample_count = function(sample_count) {
  var samp_elem = $('.sample-count').text(sample_count);
  var plural = samp_elem.siblings('.plural');
  if(sample_count == 1) {
    plural.hide();
  } else {
    plural.show();
  }
}

Interface.prototype._index_samples_by_method = function(sample_list) {
  var samples_by_method = {};
  Object.keys(sample_list).sort().forEach(function(sample_id) {
    var methods = sample_list[sample_id].methods;
    methods.forEach(function(method) {
      if(!(method in samples_by_method))
        samples_by_method[method] = {}; // Use object for O(1) lookup
      samples_by_method[method][sample_id] = true;
    });
  });
  this._samples_by_method = samples_by_method;
}

Interface.prototype._filter = function() {
  var active_methods = d3.select('#global-method-filter')
                         .selectAll('.method-choice.active')
                         .data();

  var elems = $('#sample-selector-extended tbody').find('tr');
  elems.hide();

  var self = this;
  var visible = elems.filter(function() {
    if(active_methods.length === 0)
      return false;

    var row = $(this);
    var sampid = row.find('.sampid').text().toLowerCase();
    var tumor_type = row.find('.tumor-type').text().toLowerCase();

    for(var i = 0; i < active_methods.length; i++) {
      var method = active_methods[i];
      if(!(sampid in self._samples_by_method[method]))
        return false;
    }

    var sample_filter = $('#sample-filter');
    var filter_text = $('#sample-filter').val().toLowerCase();
    if(sampid.indexOf(filter_text) === -1 && tumor_type.indexOf(filter_text) === -1)
      return false;

    return true;
  }).show();

  this._update_sample_count(visible.length);
}

Interface.prototype._activate_filters = function(sample_list) {
  var methods = Object.keys(this._samples_by_method).sort();
  var self = this;

  this.make_methods_filter('#global-method-filter', methods, function() {
    self._filter();
  });

  $('#sample-filter').keyup(function(evt) {
    self._filter();
  });
}

Interface.prototype.pick_colours = function() {
  return {
    //consensus: '#1b9e77',
    consensus: '#000000',
    broad: '#d95f02',
    mustonen095: '#7570b3',
    peifer: '#e7298a',
    theta_diploid: '#66a61e',
    vanloo_wedge: '#e6ab02'
  };
}


function main() {
  d3.json("data/index.json", function(error, sample_list) {
    if(error) return console.warn(error);
    d3.json("data/metadata.json", function(error, metadata) {
      new Interface(sample_list, metadata);
      $('#sample-selector-extended').modal('show');
    });
  });
}

main();
