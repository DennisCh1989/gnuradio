/* -*- c++ -*- */
/*
 * Copyright 2012,2014-2015 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "waterfall_sink_pr_impl.h"
#include <gnuradio/io_signature.h>
#include <gnuradio/prefs.h>
#include <string.h>
#include <volk/volk.h>
#include <iostream>

namespace gr {
  namespace qtgui {

    waterfall_sink_pr::sptr
    waterfall_sink_pr::make(const std::string &name,
			   int nconnections,
			   QWidget *parent)
    {
      return gnuradio::get_initial_sptr
	(new waterfall_sink_pr_impl(name,
				   nconnections,
				   parent));
    }

    waterfall_sink_pr_impl::waterfall_sink_pr_impl(const std::string &name,
						 int nconnections,
						 QWidget *parent)
      : sync_block("waterfall_sink_pr",
                   io_signature::make(0, nconnections, sizeof(float)),
                   io_signature::make(0, 0, 0)),
	d_name(name),
	d_nconnections(nconnections), d_nrows(200),
        d_parent(parent),
        d_port(pmt::mp("freq"))
    {
      // Required now for Qt; argc must be greater than 0 and argv
      // must have at least one valid character. Must be valid through
      // life of the qApplication:
      // http://harmattan-dev.nokia.com/docs/library/html/qt4/qapplication.html
      d_argc = 1;
      d_argv = new char;
      d_argv[0] = '\0';

      d_main_gui = NULL;

      // Perform fftshift operation;
      // this is usually desired when plotting
      
            
      d_index = 0;
      // save the last "connection" for the PDU memory
      
      d_magbufs.push_back(d_pdu_magbuf);
      memset(d_residbufs[d_nconnections], 0, d_fftsize*sizeof(float));

      initialize();

      // setup output message port to post frequency when display is
      // double-clicked
      message_port_register_out(d_port);
      message_port_register_in(d_port);
      set_msg_handler(d_port,
                      boost::bind(&waterfall_sink_f_impl::handle_set_freq, this, _1));

      // setup PDU handling input port
      message_port_register_in(pmt::mp("in"));
      set_msg_handler(pmt::mp("in"),
                      boost::bind(&waterfall_sink_f_impl::handle_pdus, this, _1));
    }

    waterfall_sink_pr_impl::~waterfall_sink_pr_impl()
    {
      if(!d_main_gui->isClosed())
        d_main_gui->close();

      for(int i = 0; i < (int)d_residbufs.size(); i++){
 	volk_free(d_magbufs[i]);
      }

      delete d_argv;
    }

    bool
    waterfall_sink_pr_impl::check_topology(int ninputs, int noutputs)
    {
      return ninputs == d_nconnections;
    }

    void
    waterfall_sink_pr_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
    {
      unsigned int ninputs = ninput_items_required.size();
      for (unsigned int i = 0; i < ninputs; i++) {
 	ninput_items_required[i] = 8191;
      }
    }

    void
    waterfall_sink_pr_impl::initialize()
    {
      if(qApp != NULL) {
 	d_qApplication = qApp;
      }
      else {
#if QT_VERSION >= 0x040500
        std::string style = prefs::singleton()->get_string("qtgui", "style", "raster");
        QApplication::setGraphicsSystem(QString(style.c_str()));
#endif
 	d_qApplication = new QApplication(d_argc, &d_argv);
      }

      // If a style sheet is set in the prefs file, enable it here.
      check_set_qss(d_qApplication);

      int numplots = (d_nconnections > 0) ? d_nconnections : 1;
      d_main_gui = new WaterfallDisplayForm(numplots, d_parent);
     
      if(d_name.size() > 0)
        set_title(d_name);

      // initialize update time to 10 times a second
      set_update_time(0.1);
    }

    void
    waterfall_sink_pr_impl::exec_()
    {
      d_qApplication->exec();
    }

    QWidget*
    waterfall_sink_pr_impl::qwidget()
    {
      return d_main_gui;
    }

#ifdef ENABLE_PYTHON
    PyObject*
    waterfall_sink_pr_impl::pyqwidget()
    {
      PyObject *w = PyLong_FromVoidPtr((void*)d_main_gui);
      PyObject *retarg = Py_BuildValue("N", w);
      return retarg;
    }
#else
    void *
    waterfall_sink_pr_impl::pyqwidget()
    {
      return NULL;
    }
#endif

    void
    waterfall_sink_pr_impl::clear_data()
    {
      d_main_gui->clearData();
    }

    
    
    void
    waterfall_sink_pr_impl::set_intensity_range(const double min,
					       const double max)
    {
      d_main_gui->setIntensityRange(min, max);
    }

    void
    waterfall_sink_pr_impl::set_update_time(double t)
    {
      //convert update time to ticks
      gr::high_res_timer_type tps = gr::high_res_timer_tps();
      d_update_time = t * tps;
      d_main_gui->setUpdateTime(t);
      d_last_time = 0;
    }

    void
    waterfall_sink_pr_impl::set_title(const std::string &title)
    {
      d_main_gui->setTitle(title.c_str());
    }

    void
    waterfall_sink_pr_impl::set_time_title(const std::string &title)
    {
        d_main_gui->setTimeTitle(title);
    }

    void
    waterfall_sink_pr_impl::set_line_label(int which, const std::string &label)
    {
      d_main_gui->setLineLabel(which, label.c_str());
    }

    void
    waterfall_sink_pr_impl::set_color_map(int which, const int color)
    {
      d_main_gui->setColorMap(which, color);
    }

    void
    waterfall_sink_pr_impl::set_line_alpha(int which, double alpha)
    {
      d_main_gui->setAlpha(which, (int)(255.0*alpha));
    }

    void
    waterfall_sink_pr_impl::set_size(int width, int height)
    {
      d_main_gui->resize(QSize(width, height));
    }

    void
    waterfall_sink_pr_impl::set_plot_pos_half(bool half)
    {
      d_main_gui->setPlotPosHalf(half);
    }

    std::string
    waterfall_sink_pr_impl::title()
    {
      return d_main_gui->title().toStdString();
    }

    std::string
    waterfall_sink_pr_impl::line_label(int which)
    {
      return d_main_gui->lineLabel(which).toStdString();
    }

    int
    waterfall_sink_pr_impl::color_map(int which)
    {
      return d_main_gui->getColorMap(which);
    }

    double
    waterfall_sink_pr_impl::line_alpha(int which)
    {
      return (double)(d_main_gui->getAlpha(which))/255.0;
    }

    void
    waterfall_sink_pr_impl::auto_scale()
    {
      d_main_gui->autoScale();
    }

    double
    waterfall_sink_pr_impl::min_intensity(int which)
    {
      return d_main_gui->getMinIntensity(which);
    }

    double
    waterfall_sink_pr_impl::max_intensity(int which)
    {
      return d_main_gui->getMaxIntensity(which);
    }

    void
    waterfall_sink_pr_impl::enable_menu(bool en)
    {
      d_main_gui->enableMenu(en);
    }

    void
    waterfall_sink_pr_impl::enable_grid(bool en)
    {
      d_main_gui->setGrid(en);
    }

    void
    waterfall_sink_pr_impl::enable_axis_labels(bool en)
    {
        d_main_gui->setAxisLabels(en);
    }

    void
    waterfall_sink_pr_impl::disable_legend()
    {
      d_main_gui->disableLegend();
    }

    
    
    void
    waterfall_sink_pr_impl::check_clicked()
    {
      if(d_main_gui->checkClicked()) {
        double freq = d_main_gui->getClickedFreq();
        message_port_pub(d_port,
                         pmt::cons(d_port,
                                   pmt::from_double(freq)));
      }
    }

    void
    waterfall_sink_pr_impl::handle_set_freq(pmt::pmt_t msg)
    {
      if(pmt::is_pair(msg)) {
        pmt::pmt_t x = pmt::cdr(msg);
        if(pmt::is_real(x)) {
          d_center_freq = pmt::to_double(x);
          d_qApplication->postEvent(d_main_gui,
                                    new SetFreqEvent(d_center_freq, d_bandwidth));
        }
      }
    }


    int
    waterfall_sink_pr_impl::work(int noutput_items,
				gr_vector_const_void_star &input_items,
				gr_vector_void_star &output_items)
    {
      int j=0;
      const float *in = (const float*)input_items[0];

      // Update the FFT size from the application
      check_clicked();

      for(int i=0; i < noutput_items; i+=d_fftsize) {
	unsigned int datasize = noutput_items - i;
	unsigned int resid = d_fftsize-d_index;

	// If we have enough input for one full FFT, do it
	if(datasize >= resid) {

	  if(gr::high_res_timer_now() - d_last_time > d_update_time) {
            for(int n = 0; n < d_nconnections; n++) {
              // Fill up residbuf with d_fftsize number of items
              in = (const float*)input_items[n];
              memcpy(d_residbufs[n]+d_index, &in[j], sizeof(float)*resid);

              for(int x = 0; x < d_fftsize; x++) {
                d_magbufs[n][x] = (double)((1.0-d_fftavg)*d_magbufs[n][x] + (d_fftavg)*d_fbuf[x]);
              }
              //volk_32f_convert_64f(d_magbufs[n], d_fbuf, d_fftsize);
            }

	    d_last_time = gr::high_res_timer_now();
	    d_qApplication->postEvent(d_main_gui,
				      new WaterfallUpdateEvent(d_magbufs,
							       d_fftsize,
							       d_last_time));
	  }

	  d_index = 0;
	  j += resid;
	}
	// Otherwise, copy what we received into the residbuf for next time
	else {
	  for(int n = 0; n < d_nconnections; n++) {
	    in = (const float*)input_items[n];
	    memcpy(d_residbufs[n]+d_index, &in[j], sizeof(float)*datasize);
	  }
	  d_index += datasize;
	  j += datasize;
	}
      }

      return j;
    }

    void
    waterfall_sink_f_impl::handle_pdus(pmt::pmt_t msg)
    {
      size_t len;
      size_t start = 0;
      pmt::pmt_t dict, samples;

      // Test to make sure this is either a PDU or a uniform vector of
      // samples. Get the samples PMT and the dictionary if it's a PDU.
      // If not, we throw an error and exit.
      if(pmt::is_pair(msg)) {
        dict = pmt::car(msg);
        samples = pmt::cdr(msg);

        pmt::pmt_t start_key = pmt::string_to_symbol("start");
        if(pmt::dict_has_key(dict, start_key)) {
          start = pmt::to_uint64(pmt::dict_ref(dict, start_key, pmt::PMT_NIL));
        }
      }
      else if(pmt::is_uniform_vector(msg)) {
        samples = msg;
      }
      else {
        throw std::runtime_error("time_sink_c: message must be either "
                                 "a PDU or a uniform vector of samples.");
      }

      len = pmt::length(samples);

      const float *in;
      if(pmt::is_f32vector(samples)) {
        in = (const float*)pmt::f32vector_elements(samples, len);
      }
      else {
        throw std::runtime_error("waterfall sink: unknown data type "
                                 "of samples; must be float.");
      }

      // Plot if we're past the last update time
      if(gr::high_res_timer_now() - d_last_time > d_update_time) {
        d_last_time = gr::high_res_timer_now();

        // Update the FFT size from the application
        fftresize();
        windowreset();
        check_clicked();

        gr::high_res_timer_type ref_start = (uint64_t)start * (double)(1.0/d_bandwidth) * 1000000;

        int stride = std::max(0, (int)(len - d_fftsize)/(int)(d_nrows));

        set_time_per_fft(1.0/d_bandwidth * stride);
        std::ostringstream title("");
        title << "Time (+" << (uint64_t)ref_start << "us)";
        set_time_title(title.str());

        int j = 0;
        size_t min = 0;
        size_t max = std::min(d_fftsize, static_cast<int>(len));
        for(size_t i=0; j < d_nrows; i+=stride) {
          // Clear residbufs if len < d_fftsize
          memset(d_residbufs[d_nconnections], 0x00, sizeof(float)*d_fftsize);

          // Copy in as much of the input samples as we can
          memcpy(d_residbufs[d_nconnections], &in[min], sizeof(float)*(max-min));

          // Apply the window and FFT; copy data into the PDU
          // magnitude buffer.
          fft(d_fbuf, d_residbufs[d_nconnections], d_fftsize);
          for(int x = 0; x < d_fftsize; x++) {
            d_pdu_magbuf[j * d_fftsize + x] = (double)d_fbuf[x];
          }

          // Increment our indices; set max up to the number of
          // samples in the input PDU.
          min += stride;
          max = std::min(max + stride, len);
          j++;
        }

        //update gui per-pdu
        d_qApplication->postEvent(d_main_gui,
                                  new WaterfallUpdateEvent(d_magbufs,
                                                           d_fftsize*d_nrows,
                                                           0));
      }
    }

  } /* namespace qtgui */
} /* namespace gr */
