/*

  An implementation of a pipe
    -- temprorary unused, being prepared for pipe processing.

  Copyright © 2018 - 2026 F.Hroch (hroch@physics.muni.cz)

  This file is part of Munipack.

  Munipack is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Munipack is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "mprocess.h"
#include <wx/wx.h>
#include <wx/app.h>
#include <wx/log.h>
#include <wx/event.h>

MuniPipe::MuniPipe(wxEvtHandler *h): parent(h), exitcode(1)
{
  Bind(wxEVT_END_PROCESS,&MuniPipe::OnFinish,this);
}

MuniPipe::~MuniPipe()
{
  wxASSERT(procs.empty());
  Unbind(wxEVT_END_PROCESS,&MuniPipe::OnFinish,this);
}

void MuniPipe::push(MuniProcess *p)
{
  if( p )
    procs.push(p);
}

void MuniPipe::Start()
{
  MuniProcess *p = procs.front();
  if( p ) {
    p->Write(OutputBuffer);
    p->OnPreProcess();
    p->OnStart();
  }
}

void MuniPipe::Stop()
{
  if( ! procs.empty() )
    procs.front()->Kill();
}

void MuniPipe::OnFinish(wxProcessEvent& event)
{
  //  wxLogDebug("MuniPipe::OnFinish %d %d",event.GetPid(),event.GetExitCode());
  exitcode = event.GetExitCode();

  if( ! procs.empty() ) {
    OutputBuffer = procs.front()->GetOutput();
    ErrorBuffer = procs.front()->GetErrors();
  }

  if( event.GetExitCode() != 0 ) {
    while( ! procs.empty() ) {
      MuniProcess *p = procs.front();
      delete p;
      procs.pop();
    }
    wxASSERT(procs.empty());
  }
  else {
    wxASSERT(!procs.empty());
    procs.front()->OnPostProcess();
    MuniProcess *p = procs.front();
    delete p;
    procs.pop();
  }

  // wxLogDebug("MuniPipe::OnFinish Remaning processes: %d",(int) procs.size());

  if( ! procs.empty() )
    Start();

  else if( parent )
    wxQueueEvent(parent,event.Clone());

  else
    event.Skip();
}

wxArrayString MuniPipe::GetOutput() const
{
  if( ! procs.empty() )
    return procs.front()->GetOutput();
  else
    return OutputBuffer;
}

wxArrayString MuniPipe::GetErrors() const
{
  if( ! procs.empty() )
    return procs.front()->GetErrors();
  else
    return ErrorBuffer;
}
